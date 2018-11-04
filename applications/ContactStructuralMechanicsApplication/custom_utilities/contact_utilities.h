// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"

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
    
/** 
 * @class ContactUtilities 
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class includes some utilities used for contact computations
 * @author Vicente Mataix Ferrandiz
 */
class ContactUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MortarUtilities
    KRATOS_CLASS_POINTER_DEFINITION( ContactUtilities );

    // Some geometrical definitions
    typedef Node<3>                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;

    /// The containers of the components of the model parts
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief This function computes the relative size of the mesh
     * @param rThisModelPart The modelpart to compute
     */

    static inline double CalculateRelativeSizeMesh(ModelPart& rThisModelPart)
    {
        return CalculateMaxNodalH(rThisModelPart)/CalculateMinimalNodalH(rThisModelPart);
    }

    /**
     * @brief This method computes the maximal nodal H
     * @param rThisModelPart The modelpart to compute
     */
    static inline double CalculateMaxNodalH(ModelPart& rThisModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& nodes_array = rThisModelPart.Nodes();

//         // Creating the max auxiliar value
//         double max_value = 0.0;
//         #pragma omp parallel for reduction(max:max_value)
//         for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
//             auto it_node = nodes_array.begin() + i;
//             KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
//             max_value = std::max(max_value, it_node->FastGetSolutionStepValue(NODAL_H));
//         }
//
//         return max_value;

        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<double> max_vector(num_threads, 0.0);
        double nodal_h;
        #pragma omp parallel for private(nodal_h)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
            nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);

            const int id = OpenMPUtils::ThisThread();

            if (nodal_h > max_vector[id])
                max_vector[id] = nodal_h;
        }

        return *std::max_element(max_vector.begin(), max_vector.end());
    }

    /**
     * @brief This method computes the mean nodal H
     * @param rThisModelPart The modelpart to compute
     */
    static inline double CalculateMeanNodalH(ModelPart& rThisModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& nodes_array = rThisModelPart.Nodes();

        double sum_nodal_h = 0.0;

        #pragma omp parallel for reduction(+:sum_nodal_h)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
            sum_nodal_h += it_node->FastGetSolutionStepValue(NODAL_H);;
        }

        return sum_nodal_h/static_cast<double>(nodes_array.size());
    }

    /**
     * @brief This method computes the minimal nodal H
     * @param rThisModelPart The modelpart to compute
     */
    static inline double CalculateMinimalNodalH(ModelPart& rThisModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& nodes_array = rThisModelPart.Nodes();

//         // Creating the min auxiliar value
//         double min_value = 0.0;
//         #pragma omp parallel for reduction(min:min_value)
//         for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
//             auto it_node = nodes_array.begin() + i;
//             KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
//             min_value = std::min(min_value, it_node->FastGetSolutionStepValue(NODAL_H));
//         }
//
//         return min_value;

        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<double> min_vector(num_threads, 0.0);
        double nodal_h;
        #pragma omp parallel for private(nodal_h)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
            nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);

            const int id = OpenMPUtils::ThisThread();

            if (nodal_h > min_vector[id])
                min_vector[id] = nodal_h;
        }

        return *std::min_element(min_vector.begin(), min_vector.end());
    }

    /**
     * @brief This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale The point to scale
     * @param Normal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
     */
    
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& PointToScale,
        const array_1d<double, 3>& Normal,
        const double LengthSearch
        )
    {        
        PointToScale.Coordinates() = PointToScale.Coordinates() + Normal * LengthSearch;
    }
    
    /**
     * @brief Calculates the distance between nodes
     * @param PointOrigin The first node
     * @param PointDestiny The second node
     */
    
    static inline double DistancePoints(
        const GeometryType::CoordinatesArrayType& PointOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny
        )
    {
        return std::sqrt((PointOrigin[0] - PointDestiny[0]) * (PointOrigin[0] - PointDestiny[0])
                       + (PointOrigin[1] - PointDestiny[1]) * (PointOrigin[1] - PointDestiny[1])
                       + (PointOrigin[2] - PointDestiny[2]) * (PointOrigin[2] - PointDestiny[2]));
    }
    
    /**
     * @brief It calculates the center updated in u_n+1 or u_n+1/2
     * @param rThisModelPart The modelpart to update
     * @param DeltaTime The increment of time considered
     * @param HalfJump If the jumpt is just half dt
     */
    
    static inline void ComputeStepJump(
        ModelPart& rThisModelPart,
        const double DeltaTime,
        const bool HalfJump = true
        )
    {
        // Time constants 
        const double velocity_constant = HalfJump ? 0.25 : 0.5;     
        const double acceleration_constant = HalfJump ? 0.125 : 0.5;
        
        // Iterate over the nodes
        NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
        // We compute the half jump
        array_1d<double, 3> new_delta_disp;
        #pragma omp parallel for private(new_delta_disp)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i)  {
            auto it_node = nodes_array.begin() + i;
            const array_1d<double, 3>& current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3>& previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            new_delta_disp = velocity_constant * DeltaTime * (current_velocity + previous_velocity) + acceleration_constant * std::pow(DeltaTime, 2) * previous_acceleration;
            if (it_node->IsFixed(DISPLACEMENT_X)) new_delta_disp[0] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Y)) new_delta_disp[1] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Z)) new_delta_disp[2] = 0.0;
            it_node->SetValue(DELTA_COORDINATES, new_delta_disp);
        }
    }
    
    /**
     * @brief It calculates the center updated in u_n+1/2
     * @param ThisGeometry The geometry to calculate
     * @return point: The center in u_n+1/2 (Newmark)
     */
    
    static inline array_1d<double, 3> GetHalfJumpCenter(GeometryType& ThisGeometry)
    {
        array_1d<double, 3> center = (ThisGeometry.Center()).Coordinates();
        
        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;
        
        // Get shape functions
        ThisGeometry.PointLocalCoordinates( local_point, center );
        ThisGeometry.ShapeFunctionsValues( N, local_point );
        
        KRATOS_DEBUG_ERROR_IF(ThisGeometry[0].Has(DELTA_COORDINATES) == false) << "WARNING:: Please call ComputeStepJump() first" << std::endl;

        const Vector new_delta_disp_center = prod(trans(GetVariableMatrix(ThisGeometry, DELTA_COORDINATES)), N);
        
        for (IndexType i = 0; i < new_delta_disp_center.size(); ++i)
            center[i] += new_delta_disp_center[i];
        
        return center;
    }
    
private:

    /**
     * @brief It calculates the matrix of a variable of a geometry
     * @param Nodes The geometry to calculate
     * @param rVarName The name of the variable to calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */

    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */
        const SizeType num_nodes = Nodes.size();
        const SizeType dim = Nodes.WorkingSpaceDimension();
        Matrix var_matrix(num_nodes, dim);

        for (IndexType i_node = 0; i_node < num_nodes; i_node++) {
            const array_1d<double, 3> value = Nodes[i_node].GetValue(rVarName);
            for (IndexType i_dof = 0; i_dof < dim; i_dof++)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

};// class ContactUtilities

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
