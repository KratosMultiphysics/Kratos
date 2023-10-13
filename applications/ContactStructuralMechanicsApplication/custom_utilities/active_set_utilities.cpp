// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_utilities/active_set_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "contact_structural_mechanics_application_variables.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
namespace ActiveSetUtilities
{
std::size_t ComputePenaltyFrictionlessActiveSet(ModelPart& rModelPart)
{
    // Defining the convergence
    IndexType is_converged = 0;

    // We get the process info
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
    if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
        const double common_epsilon = r_process_info[INITIAL_PENALTY];
        is_converged = block_for_each<SumReduction<IndexType>>(rModelPart.GetSubModelPart("Contact").Nodes(), [&](Node& rNode) {
            if (rNode.Is(SLAVE)) {
                const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;
                const double augmented_normal_pressure = epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

                rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    if (rNode.IsNot(ACTIVE)) {
                        rNode.Set(ACTIVE, true);
                        return 1;
                    }
                } else {
                    if (rNode.Is(ACTIVE)) {
                        rNode.Set(ACTIVE, false);
                        return 1;
                    }
                }
            }
            return 0;
        });
    }

    return is_converged;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<std::size_t, 2> ComputePenaltyFrictionalActiveSet(
    ModelPart& rModelPart,
    const bool PureSlip,
    const SizeType EchoLevel
    )
{
    // Auxiliary zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // Defining the convergence
    array_1d<std::size_t, 2> is_converged;
    is_converged[0] = 0;
    is_converged[1] = 0;

    const std::size_t aux_1_index = 1; 
    std::size_t& is_converged_0 = is_converged[0];
    std::size_t& is_converged_1 = is_converged[1];

    // We get the process info
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
    if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
        const double common_epsilon = r_process_info[INITIAL_PENALTY];
        const double tangent_factor = r_process_info[TANGENT_FACTOR];

        auto& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();

        block_for_each(r_nodes_array, [&](Node& rNode) {
            if (rNode.Is(SLAVE)) {
                const bool is_slip = rNode.Is(SLIP);
                const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;

                const double augmented_normal_pressure = epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

//                 // BEGIN Adding debugging value of NORMAL_GAP
//                 rNode.SetValue(NORMAL_GAP, rNode.FastGetSolutionStepValue(WEIGHTED_GAP)/rNode.GetValue(NODAL_AREA));
//                 // END Adding debugging value of NORMAL_GAP

                rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure);

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    // The friction coefficient
                    const double mu = rNode.GetValue(FRICTION_COEFFICIENT);

                    // The weighted slip
                    const array_1d<double, 3>& r_gt = rNode.FastGetSolutionStepValue(WEIGHTED_SLIP);

//                     // BEGIN Adding debugging value of TANGENT_SLIP
//                     rNode.SetValue(TANGENT_SLIP, r_gt/rNode.GetValue(NODAL_AREA));
//                     // END Adding debugging value of TANGENT_SLIP

                    // Computing the augmented tangent pressure
                    const array_1d<double,3> augmented_tangent_pressure_components = tangent_factor * epsilon * r_gt;

                    // Finally we assign and compute the norm
                    rNode.SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_tangent_pressure_components);
                    const double augmented_tangent_pressure = norm_2(augmented_tangent_pressure_components);

                    // We activate the deactivated nodes and add the contribution
                    if (rNode.IsNot(ACTIVE)) {
                        rNode.Set(ACTIVE, true);
                        AtomicAdd(is_converged_0, aux_1_index);
                    }

//                     // BEGIN Adding debugging value of TANGENTIAL_CONTACT_STRESS
//                     rNode.SetValue(TANGENTIAL_CONTACT_STRESS, augmented_tangent_pressure);
//                     // END Adding debugging value of TANGENTIAL_CONTACT_STRESS

                    // Check for the slip/stick state
                    if (augmented_tangent_pressure <= - mu * augmented_normal_pressure) { // STICK CASE // TODO: Check the <=
//                         KRATOS_WARNING_IF("ComputePenaltyFrictionalActiveSet", norm_2(r_gt) > Tolerance) << "In case of stick should be zero, if not this means that is not properly working. Node ID: " << rNode.Id() << std::endl;
//                         noalias(rNode.FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array; // NOTE: In case of stick should be zero, if not this means that is not properly working
                        KRATOS_WARNING_IF("ComputePenaltyFrictionalActiveSet", PureSlip && EchoLevel > 0) << "This node is supposed to be on STICK state. Currently working on pure slip. Node ID: " << rNode.Id() << "\tTangent pressure: " << augmented_tangent_pressure << "\tNormal x friction coeff.: " << mu * std::abs(augmented_normal_pressure)  << std::endl;
                        if (is_slip && !PureSlip) {
                            rNode.Set(SLIP, false);
                            AtomicAdd(is_converged_1, aux_1_index);
                        }
                    } else { // SLIP CASE
                        const double norm_slip = norm_2(r_gt);
                        const array_1d<double,3> tangent_direction = r_gt/norm_slip;
                        rNode.SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, mu * augmented_normal_pressure * tangent_direction);
                        if (!is_slip && !PureSlip) {
                            rNode.Set(SLIP, true);
                            AtomicAdd(is_converged_1, aux_1_index);
                        } else if (PureSlip) {
                            rNode.Set(SLIP, true);
                        }
                    }
                } else {
                    noalias(rNode.FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array;
                    if (rNode.Is(ACTIVE)) {
                        rNode.Set(ACTIVE, false);
                        rNode.Reset(SLIP);
                        AtomicAdd(is_converged_0, aux_1_index);
                    }
                }
            }
        });
    }

    return is_converged;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ComputeALMFrictionlessActiveSet(ModelPart& rModelPart)
{
    // Defining the convergence
    IndexType is_converged = 0;

    // We get the process info
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
    if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
        const double common_epsilon = r_process_info[INITIAL_PENALTY];
        const double scale_factor = r_process_info[SCALE_FACTOR];
        is_converged = block_for_each<SumReduction<IndexType>>(rModelPart.GetSubModelPart("Contact").Nodes(), [&](Node& rNode) {
            if (rNode.Is(SLAVE)) {
                const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;
                const double augmented_normal_pressure = scale_factor * rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) + epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

                rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    if (rNode.IsNot(ACTIVE)) {
                        rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = augmented_normal_pressure/scale_factor;
                        rNode.Set(ACTIVE, true);
                        return 1;
                    }
                } else {
                    if (rNode.Is(ACTIVE)) {
                        rNode.Set(ACTIVE, false);
                        return 1;
                    }
                }
            }
            return 0;
        });
    }

    return is_converged;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ComputeALMFrictionlessComponentsActiveSet(ModelPart& rModelPart)
{
    // Defining the convergence
    IndexType is_converged = 0;

    // We get the process info
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
    if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
        const double common_epsilon = r_process_info[INITIAL_PENALTY];
        const double scale_factor = r_process_info[SCALE_FACTOR];

        is_converged = block_for_each<SumReduction<IndexType>>(rModelPart.GetSubModelPart("Contact").Nodes(), [&](Node& rNode) {
            if (rNode.Is(SLAVE)) {
                const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;

                const array_1d<double,3>& r_lagrange_multiplier = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);
                const double normal_r_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);

                const double augmented_normal_pressure = scale_factor * normal_r_lagrange_multiplier + epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

                rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    if (rNode.IsNot(ACTIVE)) {
                        noalias(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = r_nodal_normal * augmented_normal_pressure/scale_factor;
                        rNode.Set(ACTIVE, true);
                        return 1;
                    }
                } else {
                    if (rNode.Is(ACTIVE)) {
                        rNode.Set(ACTIVE, false);
                        return 1;
                    }
                }
            }
            return 0;
        });
    }

    return is_converged;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<std::size_t, 2> ComputeALMFrictionalActiveSet(
    ModelPart& rModelPart,
    const bool PureSlip,
    const SizeType EchoLevel
    )
{
    // Auxiliary zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // Defining the convergence
    array_1d<std::size_t, 2> is_converged;
    is_converged[0] = 0;
    is_converged[1] = 0;

    const std::size_t aux_1_index = 1; 
    std::size_t& is_converged_0 = is_converged[0];
    std::size_t& is_converged_1 = is_converged[1];

    // We get the process info
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
    if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
        const double common_epsilon = r_process_info[INITIAL_PENALTY];
        const double scale_factor = r_process_info[SCALE_FACTOR];
        const double tangent_factor = r_process_info[TANGENT_FACTOR];

        // Slip convergence enhancers
        const double slip_threshold = r_process_info.Has(SLIP_THRESHOLD) ? r_process_info[SLIP_THRESHOLD] : 0.0;
        const double slip_augmentation_coefficient = r_process_info.Has(SLIP_AUGMENTATION_COEFFICIENT) ? r_process_info[SLIP_AUGMENTATION_COEFFICIENT] : 0.0;

        auto& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();

        block_for_each(r_nodes_array, [&](Node& rNode) {
            if (rNode.Is(SLAVE)) {
                const bool is_slip = rNode.Is(SLIP);
                const double epsilon = rNode.Has(INITIAL_PENALTY) ? rNode.GetValue(INITIAL_PENALTY) : common_epsilon;

                const array_1d<double,3>& r_lagrange_multiplier = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                const array_1d<double,3>& r_nodal_normal = rNode.FastGetSolutionStepValue(NORMAL);
                const double normal_lagrange_multiplier = inner_prod(r_nodal_normal, r_lagrange_multiplier);

                const double augmented_normal_pressure = scale_factor * normal_lagrange_multiplier + epsilon * rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

                rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure);

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    // The friction coefficient
                    const double mu = rNode.GetValue(FRICTION_COEFFICIENT);

                    // The weighted slip
                    const array_1d<double, 3>& r_gt = rNode.FastGetSolutionStepValue(WEIGHTED_SLIP);

                    // Computing the augmented tangent pressure
                    const array_1d<double, 3> tangent_lagrange_multiplier = r_lagrange_multiplier - normal_lagrange_multiplier * r_nodal_normal;
                    const array_1d<double, 3> augmented_tangent_pressure_components = is_slip ? scale_factor * tangent_lagrange_multiplier + slip_augmentation_coefficient * tangent_factor * epsilon * r_gt : scale_factor * tangent_lagrange_multiplier + tangent_factor * epsilon * r_gt;

                    // Finally we assign and compute the norm
                    rNode.SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_tangent_pressure_components);
                    const double augmented_tangent_pressure = norm_2(augmented_tangent_pressure_components);

                    // We activate the deactivated nodes and add the contribution
                    if (rNode.IsNot(ACTIVE)) {
                        noalias(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = r_nodal_normal * augmented_normal_pressure/scale_factor + (mu < std::numeric_limits<double>::epsilon() ? zero_array : augmented_tangent_pressure_components/scale_factor);
                        rNode.Set(ACTIVE, true);
                        AtomicAdd(is_converged_0, aux_1_index);
                    }

                    // Check for the slip/stick state
                    const double threshold_value = is_slip ? 1.0 - slip_threshold : 1.0;
                    const bool slip_check = augmented_tangent_pressure/(- mu * augmented_normal_pressure) > threshold_value ? true : false;
                    if (!slip_check) { // STICK CASE
//                     if (augmented_tangent_pressure <= - mu * augmented_normal_pressure) { // STICK CASE // TODO: Check the <=
//                             KRATOS_WARNING_IF("ComputeALMFrictionalActiveSet", norm_2(r_gt) > Tolerance) << "In case of stick should be zero, if not this means that is not properly working. Node ID: " << rNode.Id() << std::endl;
//                             noalias(rNode.FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array; // NOTE: In case of stick should be zero, if not this means that is not properly working
                        KRATOS_WARNING_IF("ComputeALMFrictionalActiveSet", PureSlip && EchoLevel > 0) << "This node is supposed to be on STICK state. Currently working on pure slip. Node ID: " << rNode.Id() << "\tTangent pressure: " << augmented_tangent_pressure << "\tNormal x friction coeff.: " << mu * std::abs(augmented_normal_pressure)  << std::endl;
                        if (is_slip && !PureSlip) {
                            rNode.Set(SLIP, false);
                            AtomicAdd(is_converged_1, aux_1_index);
                        }
                    } else { // SLIP CASE
                        const array_1d<double, 3> tangent_direction = tangent_lagrange_multiplier/norm_2(tangent_lagrange_multiplier);
                        const array_1d<double, 3> augmented_contact_tangent_pressure = - mu * augmented_normal_pressure * tangent_direction;
                        rNode.SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, augmented_contact_tangent_pressure);
//                         noalias(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = normal_lagrange_multiplier * r_nodal_normal + augmented_contact_tangent_pressure/scale_factor;

                        if (!is_slip && !PureSlip) {
                            rNode.Set(SLIP, true);
                            AtomicAdd(is_converged_1, aux_1_index);
                        } else if (PureSlip) {
                            rNode.Set(SLIP, true);
                        }
                    }
                } else {
                    noalias(rNode.FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array;
                    if (rNode.Is(ACTIVE)) {
                        rNode.Set(ACTIVE, false);
                        rNode.Reset(SLIP);
                        AtomicAdd(is_converged_0, aux_1_index);
                    }
                }
            }
        });
    }

    return is_converged;
}

} // namespace ActiveSetUtilities
} // namespace Kratos
