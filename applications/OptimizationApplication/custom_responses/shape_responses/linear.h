// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef LINEAR_RESPONSE_H
#define LINEAR_RESPONSE_H

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

class KRATOS_API(OPTIMIZATION_APPLICATION) Linear : public Response
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

    /// Pointer definition of Linear
    KRATOS_CLASS_POINTER_DEFINITION(Linear);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Linear(std::string ResponseName, Model& rModel, Parameters ResponseSettings )
        : Response(ResponseName,"shape",ResponseSettings),mrModel(rModel){}

    /// Destructor.
    virtual ~Linear()
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

        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            // check if control_obj has surface condition and root_model_part has elements
            KRATOS_ERROR_IF_NOT(r_eval_object.Conditions().size()>0)
            <<"Linear::Initialize: evaluated object "<<eval_obj.GetString()<<" must have surface conditions !"<<std::endl;
        }

    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        return 1.0;
    };    
    // --------------------------------------------------------------------------
    void CalculateGradient() override {

        for(auto& eval_obj : mResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            VariableUtils().SetHistoricalVariableToZero(D_LINEAR_D_X, r_eval_object.Nodes());
            for(auto& cond_i : r_eval_object.Conditions()){
                auto& r_geom = cond_i.GetGeometry();
                unsigned int dimension = r_geom.WorkingSpaceDimension();
                unsigned int number_of_nodes = r_geom.size();
                unsigned int mat_size = dimension * number_of_nodes;
                Matrix MassMatrix = ZeroMatrix( mat_size, mat_size );


                const array_1d<double,3> local_coords = ZeroVector(3);
                array_1d<double,3> cond_normal = r_geom.UnitNormal(local_coords);
                const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
                const auto& integration_points = r_geom.IntegrationPoints( integration_method );
                const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method); 
                Matrix J0(dimension, dimension);


                array_1d<double, 3 > tangent_xi, tangent_eta;
                Matrix J(3, 2);
                for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

                    r_geom.Jacobian(J, point_number, integration_method);
                    const double detJ = MathUtils<double>::GeneralizedDet(J);
                    const double integration_weight = integration_points[point_number].Weight() * detJ;            
                    const auto& rN = row(Ncontainer, point_number);


                    tangent_xi[0]  = J(0, 0);
                    tangent_eta[0] = J(0, 1);
                    tangent_xi[1]  = J(1, 0);
                    tangent_eta[1] = J(1, 1);
                    tangent_xi[2]  = J(2, 0);
                    tangent_eta[2] = J(2, 1);

                    array_1d<double, 3 > normal;
                    MathUtils<double>::UnitCrossProduct(normal, tangent_eta, tangent_xi);
                    // Calculating the pressure on the gauss point
                    double pressure = 0.0;
                    for (std::size_t ii = 0; ii < number_of_nodes; ii++) {
                        pressure += rN[ii] * 1.0;
                    }

                    normal(0) = 0.0;
                    normal(1) = 1.0;
                    normal(2) = 0.0;

                    int node_index = 0;
                    for (auto& node_i : cond_i.GetGeometry()){
                        array_3d gradient(3, 0.0);
                        gradient = pressure * rN[node_index] * integration_weight * normal;
                        noalias(node_i.FastGetSolutionStepValue(D_LINEAR_D_X)) += gradient;
                        node_index++;
                    }
                }
            } 
            // for(auto& node_i : r_eval_object.Nodes()){

            //     array_3d gradient(3, 0.0);
            //     gradient(1) = 1;
            //     noalias(node_i.FastGetSolutionStepValue(D_LINEAR_D_X)) = gradient;
            // }          
        }
    };  

    // --------------------------------------------------------------------------
       
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
        return "Linear";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Linear";
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
//      Linear& operator=(Linear const& rOther);

    /// Copy constructor.
//      Linear(Linear const& rOther);


    ///@}

}; // Class Linear

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // LINEAR_RESPONSE_H
