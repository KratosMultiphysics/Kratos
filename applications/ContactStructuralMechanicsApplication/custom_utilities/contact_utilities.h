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

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
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
    typedef Point::CoordinatesArrayType              CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;

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
     * @param rModelPart The modelpart to compute
     */
    static inline double CalculateRelativeSizeMesh(ModelPart& rModelPart)
    {
        return CalculateMaxNodalH(rModelPart)/CalculateMinimalNodalH(rModelPart);
    }

    /**
     * @brief This method computes the maximal nodal H
     * @param rModelPart The modelpart to compute
     */
    static inline double CalculateMaxNodalH(ModelPart& rModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();

//         // Creating the max auxiliar value
//         double max_value = 0.0;
//
//         #pragma omp parallel for reduction(max:max_value)
//         for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
//             auto it_node = it_node_begin + i;
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
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
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
     * @param rModelPart The modelpart to compute
     */
    static inline double CalculateMeanNodalH(ModelPart& rModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        // Creating the sum auxiliar value
        double sum_nodal_h = 0.0;

        #pragma omp parallel for reduction(+:sum_nodal_h)
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not added" << std::endl;
            sum_nodal_h += it_node->FastGetSolutionStepValue(NODAL_H);;
        }

        return sum_nodal_h/static_cast<double>(r_nodes_array.size());
    }

    /**
     * @brief This method computes the minimal nodal H
     * @param rModelPart The modelpart to compute
     */
    static inline double CalculateMinimalNodalH(ModelPart& rModelPart)
    {
        // We iterate over the nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();

//         // Creating the min auxiliar value
//         double min_value = 0.0;
//
//         #pragma omp parallel for reduction(min:min_value)
//         for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
//             auto it_node = it_node_begin + i;
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
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
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
     * @param rPointToScale The point to scale
     * @param rNormal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
     */
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& rPointToScale,
        const array_1d<double, 3>& rNormal,
        const double LengthSearch
        )
    {
        noalias(rPointToScale.Coordinates()) = rPointToScale.Coordinates() + rNormal * LengthSearch;
    }

    /**
     * @brief Calculates the distance between nodes
     * @param rPointOrigin The first node
     * @param rPointDestiny The second node
     */
    static inline double DistancePoints(
        const GeometryType::CoordinatesArrayType& rPointOrigin,
        const GeometryType::CoordinatesArrayType& rPointDestiny
        )
    {
        return std::sqrt((rPointOrigin[0] - rPointDestiny[0]) * (rPointOrigin[0] - rPointDestiny[0])
                       + (rPointOrigin[1] - rPointDestiny[1]) * (rPointOrigin[1] - rPointDestiny[1])
                       + (rPointOrigin[2] - rPointDestiny[2]) * (rPointOrigin[2] - rPointDestiny[2]));
    }

    /**
     * @brief It calculates the center updated in u_n+1 or u_n+1/2
     * @param rModelPart The modelpart to update
     * @param DeltaTime The increment of time considered
     * @param HalfJump If the jumpt is just half dt
     */
    static inline void ComputeStepJump(
        ModelPart& rModelPart,
        const double DeltaTime,
        const bool HalfJump = true
        )
    {
        // Time constants
        const double velocity_constant = HalfJump ? 0.25 : 0.5;
        const double acceleration_constant = HalfJump ? 0.125 : 0.5;

        // Iterate over the nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();

        // Node iterator
        const auto it_node_begin = r_nodes_array.begin();

        // We compute the half jump
        array_1d<double, 3> new_delta_disp;
        #pragma omp parallel for firstprivate(new_delta_disp)
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i)  {
            auto it_node = it_node_begin + i;
            const array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3>& r_previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            noalias(new_delta_disp) = velocity_constant * DeltaTime * (r_current_velocity + r_previous_velocity) + acceleration_constant * std::pow(DeltaTime, 2) * r_previous_acceleration;
            if (it_node->IsFixed(DISPLACEMENT_X)) new_delta_disp[0] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Y)) new_delta_disp[1] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Z)) new_delta_disp[2] = 0.0;
            it_node->SetValue(DELTA_COORDINATES, new_delta_disp);
        }
    }

    /**
     * @brief It checks the activity of the current contact simulation
     * @param rModelPart The modelpart to check the activity
     * @param ThrowError If an error is thrown
     */
    static inline bool CheckActivity(
        ModelPart& rModelPart,
        const bool ThrowError = true
        )
    {
        // Iterate over the nodes
        NodesArrayType& r_nodes_array = rModelPart.Nodes();

        // Node iterator
        const auto it_node_begin = r_nodes_array.begin();

        // We compute the half jump
        IndexType aux_check = 0;
        #pragma omp parallel for reduction(+:aux_check)
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i)  {
            auto it_node = it_node_begin + i;
            if (it_node->Is(SLAVE)) {
                if (it_node->Is(ACTIVE)) {
                    aux_check += 1;
                }
            }
        }

        const bool is_active = aux_check == 0 ?  false : true;

        KRATOS_ERROR_IF(ThrowError && !is_active) << "CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;

        return is_active;
    }

    /**
     * @brief It computes the explicit contributions of the conditions
     * @param rModelPart The modelpart to update
     */
    static inline void ComputeExplicitContributionConditions(ModelPart& rModelPart)
    {
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        KRATOS_TRACE_IF("Empty model part", r_conditions_array.size() == 0) << "YOUR COMPUTING CONTACT MODEL PART IS EMPTY" << std::endl;
        const auto it_cond_begin = r_conditions_array.begin();
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->AddExplicitContribution(r_process_info);
        }
    }

    /**
     * @brief It activates the conditions with active nodes
     * @param rModelPart The modelpart to check
     */
    static inline void ActivateConditionWithActiveNodes(ModelPart& rModelPart)
    {
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        KRATOS_TRACE_IF("Empty model part", r_conditions_array.size() == 0) << "YOUR COMPUTING CONTACT MODEL PART IS EMPTY" << std::endl;
        const auto it_cond_begin = r_conditions_array.begin();

        bool is_active = false;
        #pragma omp parallel for firstprivate(is_active)
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            GeometryType& r_geometry = it_cond->GetGeometry();
            is_active = false;
            for ( IndexType i_node = 0; i_node < r_geometry.size(); ++i_node ) {
                if (r_geometry[i_node].Is(ACTIVE)) {
                    is_active = true;
                    break;
                }
            }
            it_cond->Set(ACTIVE, is_active);
        }
    }

    /**
     * @brief It calculates the center updated in u_n+1/2
     * @param rThisGeometry The geometry to calculate
     * @return point: The center in u_n+1/2 (Newmark)
     */
    static inline array_1d<double, 3> GetHalfJumpCenter(GeometryType& rThisGeometry)
    {
        array_1d<double, 3> center = (rThisGeometry.Center()).Coordinates();

        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;

        // Get shape functions
        rThisGeometry.PointLocalCoordinates( local_point, center );
        rThisGeometry.ShapeFunctionsValues( N, local_point );

        KRATOS_DEBUG_ERROR_IF_NOT(rThisGeometry[0].Has(DELTA_COORDINATES)) << "Please call ComputeStepJump() first" << std::endl;

        const Vector new_delta_disp_center = prod(trans(GetVariableMatrix(rThisGeometry, DELTA_COORDINATES)), N);

        for (IndexType i = 0; i < new_delta_disp_center.size(); ++i)
            center[i] += new_delta_disp_center[i];

        return center;
    }

private:

    /**
     * @brief It calculates the matrix of a variable of a geometry
     * @param rNodes The geometry to calculate
     * @param rVarName The name of the variable to calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
    static inline Matrix GetVariableMatrix(
        const GeometryType& rNodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */
        const SizeType num_nodes = rNodes.size();
        const SizeType dim = rNodes.WorkingSpaceDimension();
        Matrix var_matrix(num_nodes, dim);

        for (IndexType i_node = 0; i_node < num_nodes; i_node++) {
            const array_1d<double, 3> value = rNodes[i_node].GetValue(rVarName);
            for (IndexType i_dof = 0; i_dof < dim; i_dof++)
                var_matrix(i_node, i_dof) = value[i_dof];
        }

        return var_matrix;
    }

};// class ContactUtilities

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
