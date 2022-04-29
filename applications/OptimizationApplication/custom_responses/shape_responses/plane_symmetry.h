// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef PLANE_SYMMETRY_RESPONSE_H
#define PLANE_SYMMETRY_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_responses/response.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "processes/find_conditions_neighbours_process.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(OPTIMIZATION_APPLICATION) PlaneSymmetry : public Response
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef array_1d<double,3> array_3d;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of PlaneSymmetry
    KRATOS_CLASS_POINTER_DEFINITION(PlaneSymmetry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlaneSymmetry(std::string ResponseName, Model& rModel, Parameters ResponseSettings )
        : Response(ResponseName,"shape",ResponseSettings),mrModel(rModel){}

    /// Destructor.
    virtual ~PlaneSymmetry()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override {
        int total_num_points = 0;
        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            // check if control_obj has surface condition and root_model_part has elements
            KRATOS_ERROR_IF_NOT(r_eval_object.Conditions().size()>0)
            <<"PlaneSymmetry::Initialize: evaluated object "<<eval_obj.GetString()<<" must have surface conditions !"<<std::endl;
            FindConditionsNeighboursProcess find_conditions_neighbours_process(r_eval_object, r_eval_object.GetProcessInfo()[DOMAIN_SIZE]);
            find_conditions_neighbours_process.Execute(); 
            total_num_points += r_eval_object.Nodes().size();
        }

        mListOfNodesInModelPart.resize(total_num_points);
        int counter = 0;
        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());            
            for (ModelPart::NodesContainerType::iterator node_it = r_eval_object.NodesBegin(); node_it != r_eval_object.NodesEnd(); ++node_it)
            {
                NodeTypePointer pnode = *(node_it.base());
                mListOfNodesInModelPart[counter++] = pnode;
            }             
        }

        plane_normal(0) = 1.0;
        plane_normal(1) = 0.0;
        plane_normal(2) = 0.0;

        plane_point(0) = 0.0;
        plane_point(1) = 10.0;
        plane_point(2) = 5.0;
        
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {

        const size_t bucket_size = 100;
        KDTree search_tree(mListOfNodesInModelPart.begin(), mListOfNodesInModelPart.end(), bucket_size);

        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());            
            for (auto& r_node : r_eval_object.Nodes()){
                const array_1d<double, 3>& coords = r_node.Coordinates();
                array_1d<double, 3> reflected_point = 2.0 * inner_prod(plane_point-coords,plane_normal) * plane_normal + coords;
                NodeType reflected_node(0,reflected_point[0],reflected_point[1],reflected_point[2]);                
                double distance;
                NodeTypePointer p_neighbor = search_tree.SearchNearestPoint(reflected_node, distance);

                array_1d<double, 3> closest_projection_point = p_neighbor->Coordinates();
                double smallest_dist = distance;
                const GlobalPointersVector<Condition>& rConditions = p_neighbor->GetValue(NEIGHBOUR_CONDITIONS);

                for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
                {
                    // Get geometry of current condition
                    Condition rCondition = rConditions[c_itr];
                    Condition::GeometryType& geom_i = rCondition.GetGeometry();
                    const array_1d<double,3> local_coords = ZeroVector(3);
                    array_1d<double,3> cond_normal = geom_i.UnitNormal(local_coords);
                    array_1d<double, 3> a_point_on_cond = geom_i.Center();
                    array_1d<double, 3> projection_on_cond = inner_prod(a_point_on_cond-reflected_point,cond_normal) * cond_normal + reflected_point;
                    array_1d<double, 3> dist_vec = projection_on_cond - reflected_point;
                    double dist = sqrt(inner_prod(dist_vec,dist_vec));
                    if(dist<smallest_dist){
                        closest_projection_point =  projection_on_cond; 
                        smallest_dist = dist;
                    }
                      
                }
                auto& node_i_nn_reflected_point = r_node.FastGetSolutionStepValue(NEAREST_NEIGHBOUR_POINT);
                auto& node_i_nn_reflected_dist = r_node.FastGetSolutionStepValue(NEAREST_NEIGHBOUR_DIST);
                node_i_nn_reflected_point = closest_projection_point;
                node_i_nn_reflected_dist = smallest_dist;
            }
        }        

        double total_value = 0.0;
        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            VariableUtils().SetHistoricalVariableToZero(D_PLANE_SYMMETRY_D_X, r_eval_object.Nodes());
            for(auto& cond_i : r_eval_object.Conditions()){
                const auto& r_geom = cond_i.GetGeometry();
                const double g_i = CalculateConditionValue(cond_i);
                total_value += g_i;

                // Compute sensitivities using finite differencing in the three spatial direction
                array_3d gradient(3, 0.0);

                for (auto& node_i : cond_i.GetGeometry()){

                    // Apply pertubation in X-direction
                    double g_i_after_fd = 0.0;
                    node_i.X() += mDelta;
                    node_i.X0() += mDelta;
                    g_i_after_fd = CalculateConditionValue(cond_i);
                    gradient[0] = (g_i_after_fd - g_i) / mDelta;
                    node_i.X() -= mDelta;
                    node_i.X0() -= mDelta;

                    // Apply pertubation in Y-direction
                    g_i_after_fd = 0.0;
                    node_i.Y() += mDelta;
                    node_i.Y0() += mDelta;
                    g_i_after_fd = CalculateConditionValue(cond_i);
                    gradient[1] = (g_i_after_fd - g_i) / mDelta;
                    node_i.Y() -= mDelta;
                    node_i.Y0() -= mDelta;

                    // Apply pertubation in Z-direction
                    g_i_after_fd = 0.0;
                    node_i.Z() += mDelta;
                    node_i.Z0() += mDelta;
                    g_i_after_fd = CalculateConditionValue(cond_i);
                    gradient[2] = (g_i_after_fd - g_i) / mDelta;
                    node_i.Z() -= mDelta;
                    node_i.Z0() -= mDelta;

                    // Add to aggregated sensitivities
                    noalias(node_i.FastGetSolutionStepValue(D_PLANE_SYMMETRY_D_X)) += gradient;
                }


            }
                      
        }  
        return total_value;
    };    
    // --------------------------------------------------------------------------
    void CalculateGradient() override {};  

    // --------------------------------------------------------------------------
    double CalculateConditionValue(const Condition& rFace) {

        const auto& r_geom = rFace.GetGeometry();
        std::size_t dimension = r_geom.WorkingSpaceDimension();

        Vector nodal_dist_vec = ZeroVector(r_geom.PointsNumber());

        int counter = 0;
        for(auto& node_i : r_geom){

            const array_1d<double, 3>& coords = node_i.Coordinates();
            array_1d<double, 3> reflected_point = 2.0 * inner_prod(plane_point-coords,plane_normal) * plane_normal + coords;
            auto& node_i_nn_reflected_point = node_i.FastGetSolutionStepValue(NEAREST_NEIGHBOUR_POINT);
            array_1d<double, 3> dist_vec = node_i_nn_reflected_point - reflected_point;
            double dist = sqrt(inner_prod(dist_vec,dist_vec));
            nodal_dist_vec[counter] = dist;
            counter++;
        }

        Matrix J0(dimension, dimension);

        const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
        const auto& integration_points = r_geom.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        double g_i_val = 0.0;
        for ( std::size_t point_number = 0; point_number < integration_points.size(); ++point_number ) {
            GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::Det(J0);
            const double integration_weight = integration_points[point_number].Weight() * detJ0;
            const Vector& rN = row(Ncontainer,point_number);
            double value_at_gp = inner_prod(nodal_dist_vec,rN);        
            g_i_val += value_at_gp * integration_weight;
        }
    
        return g_i_val;
    };        

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "PlaneSymmetry";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PlaneSymmetry";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    // Initialized by class constructor
    Model& mrModel;
    KDTree::Pointer mpSearchTree;
    NodeVector mListOfNodesInModelPart;
    double mDelta = 1e-6;

    array_1d<double,3> plane_normal;
    array_1d<double,3> plane_point;
    
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      PlaneSymmetry& operator=(PlaneSymmetry const& rOther);

    /// Copy constructor.
//      PlaneSymmetry(PlaneSymmetry const& rOther);


    ///@}

}; // Class PlaneSymmetry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // PLANE_SYMMETRY_RESPONSE_H
