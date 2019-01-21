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

#if !defined(KRATOS_ACTIVE_SET_UTILITIES)
#define KRATOS_ACTIVE_SET_UTILITIES

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
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
 * @class ActiveSetUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class includes some utilities used for contact active set computations
 * @author Vicente Mataix Ferrandiz
 */
class ActiveSetUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MortarUtilities
    KRATOS_CLASS_POINTER_DEFINITION( ActiveSetUtilities );

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
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    static inline std::size_t ComputePenaltyFrictionlessActiveSet(ModelPart& rModelPart)
    {
        // Defining the convergence
        IndexType is_converged = 0;

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];

            NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            #pragma omp parallel for reduction(+:is_converged)
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(SLAVE)) {
                    const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;
                    const double augmented_normal_pressure = epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                    it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                    if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                        if (it_node->IsNot(ACTIVE)) {
                            it_node->Set(ACTIVE, true);
                            is_converged += 1;
                        }
                    } else {
                        if (it_node->Is(ACTIVE)) {
                            it_node->Set(ACTIVE, false);
                            is_converged += 1;
                        }
                    }
                }
            }
        }

        return is_converged;
    }

    /**
     * @brief This function computes the active set for penalty frictional cases
     * @param rThisModelPart The modelpart to compute
     */
    static inline array_1d<std::size_t, 2> ComputePenaltyFrictionalActiveSet(ModelPart& rModelPart)
    {
        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // Defining the convergence
        array_1d<std::size_t, 2> is_converged;
        is_converged[0] = 0;
        is_converged[1] = 0;

        std::size_t& is_converged_0 = is_converged[0];
        std::size_t& is_converged_1 = is_converged[1];

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];
            const double tangent_factor = r_process_info[TANGENT_FACTOR];

            NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(SLAVE)) {
                    const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;

                    const double augmented_normal_pressure = epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                    it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure);

                    if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                        // The friction coefficient
                        const double mu = it_node->GetValue(FRICTION_COEFFICIENT);

                        // The weighted slip
                        const array_1d<double, 3>& r_gt = it_node->FastGetSolutionStepValue(WEIGHTED_SLIP);

                        // Computing the augmented tangent pressure
                        const array_1d<double,3> augmented_tangent_pressure_components = tangent_factor * epsilon * r_gt;

                        // Finally we assign and compute the norm
                        it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_tangent_pressure_components);
                        const double augmented_tangent_pressure = norm_2(augmented_tangent_pressure_components);

                        // We activate the deactivated nodes and add the contribution
                        if (it_node->IsNot(ACTIVE)) {
                            it_node->Set(ACTIVE, true);
                            #pragma omp atomic
                            is_converged_0 += 1;
                        }

                        // Check for the slip/stick state
                        if (augmented_tangent_pressure <= - mu * augmented_normal_pressure) { // STICK CASE
//                             KRATOS_WARNING_IF("PenaltyFrictionalMortarConvergenceCriteria", norm_2(r_gt) > Tolerance) << "In case of stick should be zero, if not this means that is not properly working. Node ID: " << it_node->Id() << std::endl;
//                             noalias(it_node->FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array; // NOTE: In case of stick should be zero, if not this means that is not properly working
                            if (it_node->Is(SLIP)) {
                                it_node->Set(SLIP, false);
                                #pragma omp atomic
                                is_converged_1 += 1;
                            }
                        } else { // SLIP CASE
                            if (it_node->IsNot(SLIP)) {
                                it_node->Set(SLIP, true);
                                #pragma omp atomic
                                is_converged_1 += 1;
                            }
                        }
                    } else {
                        noalias(it_node->FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array;
                        if (it_node->Is(ACTIVE)) {
                            it_node->Set(ACTIVE, false);
                            it_node->Set(SLIP, false);
                            is_converged_0 += 1;
                        }
                    }
                }
            }
        }

        return is_converged;
    }

    /**
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    static inline std::size_t ComputeALMFrictionlessActiveSet(ModelPart& rModelPart)
    {
        // Defining the convergence
        IndexType is_converged = 0;

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];
            const double scale_factor = r_process_info[SCALE_FACTOR];

            NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            #pragma omp parallel for reduction(+:is_converged)
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(SLAVE)) {
                    const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;
                    const double augmented_normal_pressure = scale_factor * it_node->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) + epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                    it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                    if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                        if (it_node->IsNot(ACTIVE)) {
                            it_node->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = augmented_normal_pressure/scale_factor;
                            it_node->Set(ACTIVE, true);
                            is_converged += 1;
                        }
                    } else {
                        if (it_node->Is(ACTIVE)) {
                            it_node->Set(ACTIVE, false);
                            is_converged += 1;
                        }
                    }
                }
            }
        }

        return is_converged;
    }

    /**
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    static inline std::size_t ComputeALMFrictionlessComponentsActiveSet(ModelPart& rModelPart)
    {
        // Defining the convergence
        IndexType is_converged = 0;

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];
            const double scale_factor = r_process_info[SCALE_FACTOR];

            NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            #pragma omp parallel for reduction(+:is_converged)
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(SLAVE)) {
                    const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;

                    const array_1d<double,3>& r_lagrange_multiplier = it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    const array_1d<double,3>& r_nodal_normal = it_node->FastGetSolutionStepValue(NORMAL);
                    const double normal_r_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);

                    const double augmented_normal_pressure = scale_factor * normal_r_lagrange_multiplier + epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                    it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                    if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                        if (it_node->IsNot(ACTIVE)) {
                            noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = r_nodal_normal * augmented_normal_pressure/scale_factor;
                            it_node->Set(ACTIVE, true);
                            is_converged += 1;
                        }
                    } else {
                        if (it_node->Is(ACTIVE)) {
                            it_node->Set(ACTIVE, false);
                            is_converged += 1;
                        }
                    }
                }
            }
        }

        return is_converged;
    }

    /**
     * @brief This function computes the active set for penalty frictional cases
     * @param rThisModelPart The modelpart to compute
     */
    static inline array_1d<std::size_t, 2> ComputeALMFrictionalActiveSet(ModelPart& rModelPart)
    {
        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // Defining the convergence
        array_1d<std::size_t, 2> is_converged;
        is_converged[0] = 0;
        is_converged[1] = 0;

        std::size_t& is_converged_0 = is_converged[0];
        std::size_t& is_converged_1 = is_converged[1];

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];
            const double scale_factor = r_process_info[SCALE_FACTOR];
            const double tangent_factor = r_process_info[TANGENT_FACTOR];

            NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
            const auto it_node_begin = r_nodes_array.begin();

            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(SLAVE)) {
                    const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;

                    const array_1d<double,3>& r_lagrange_multiplier = it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    const array_1d<double,3>& r_nodal_normal = it_node->FastGetSolutionStepValue(NORMAL);
                    const double normal_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);

                    const double augmented_normal_pressure = scale_factor * normal_lagrange_multiplier + epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                    it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure);

                    if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                        // The friction coefficient
                        const double mu = it_node->GetValue(FRICTION_COEFFICIENT);

                        // The weighted slip
                        const array_1d<double, 3>& r_gt = it_node->FastGetSolutionStepValue(WEIGHTED_SLIP);

                        // Computing the augmented tangent pressure
                        const array_1d<double,3> tangent_lagrange_multiplier = r_lagrange_multiplier - normal_lagrange_multiplier * r_nodal_normal;
                        const array_1d<double,3> augmented_tangent_pressure_components = scale_factor * tangent_lagrange_multiplier + tangent_factor * epsilon * r_gt;

                        // Finally we assign and compute the norm
                        it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_tangent_pressure_components);
                        const double augmented_tangent_pressure = norm_2(augmented_tangent_pressure_components);

                        // We activate the deactivated nodes and add the contribution
                        if (it_node->IsNot(ACTIVE)) {
                            noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = r_nodal_normal * augmented_normal_pressure/scale_factor + augmented_tangent_pressure_components/scale_factor;
                            it_node->Set(ACTIVE, true);
                            #pragma omp atomic
                            is_converged_0 += 1;
                        }

                        // Check for the slip/stick state
                        if (augmented_tangent_pressure <= - mu * augmented_normal_pressure) { // STICK CASE
//                             KRATOS_WARNING_IF("PenaltyFrictionalMortarConvergenceCriteria", norm_2(r_gt) > Tolerance) << "In case of stick should be zero, if not this means that is not properly working. Node ID: " << it_node->Id() << std::endl;
//                             noalias(it_node->FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array; // NOTE: In case of stick should be zero, if not this means that is not properly working
                            if (it_node->Is(SLIP)) {
                                it_node->Set(SLIP, false);
                                #pragma omp atomic
                                is_converged_1 += 1;
                            }
                        } else { // SLIP CASE
                            if (it_node->IsNot(SLIP)) {
                                it_node->Set(SLIP, true);
                                #pragma omp atomic
                                is_converged_1 += 1;
                            }
                        }
                    } else {
                        noalias(it_node->FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array;
                        if (it_node->Is(ACTIVE)) {
                            it_node->Set(ACTIVE, false);
                            it_node->Set(SLIP, false);
                            #pragma omp atomic
                            is_converged_0 += 1;
                        }
                    }
                }
            }
        }

        return is_converged;
    }

};// class ActiveSetUtilities

}
#endif /* KRATOS_ACTIVE_SET_UTILITIES defined */
