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
#include "utilities/variable_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/compute_dynamic_factor_process.h"

namespace Kratos
{
void ComputeDynamicFactorProcess::Execute()
{
    KRATOS_TRY;

    // Getting process info
    ProcessInfo& r_process_info = mrThisModelPart.GetProcessInfo();

    // Impact time duration
    const double max_gap_factor = r_process_info.Has(MAX_GAP_FACTOR) ? r_process_info[MAX_GAP_FACTOR] : 1.0;
    const double max_gap_threshold = r_process_info.Has(MAX_GAP_THRESHOLD) ? r_process_info[MAX_GAP_THRESHOLD] : 0.0;
    const double common_epsilon = r_process_info[INITIAL_PENALTY];

    // We iterate over the node
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    block_for_each(r_nodes_array, [&max_gap_threshold, &common_epsilon, &max_gap_factor](Node& rNode) {
        // Computing only on SLAVE nodes
        if (rNode.Is(SLAVE) && rNode.Is(ACTIVE)) {
            // Weighted values
            const double normal_area  = rNode.GetValue(NODAL_AREA);
            const double current_gap  = rNode.FastGetSolutionStepValue(WEIGHTED_GAP)/normal_area;
            rNode.SetValue(NORMAL_GAP, current_gap);
            const double previous_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP, 1)/normal_area;

            // Computing actual logistic factor
            if (max_gap_threshold > 0.0 && current_gap <= 0.0) {
                const double logistic_factor = ComputeLogisticFactor(max_gap_threshold, current_gap);
                rNode.SetValue(INITIAL_PENALTY, common_epsilon * (1.0 + logistic_factor * max_gap_factor));
            } else {
                rNode.SetValue(INITIAL_PENALTY, common_epsilon);
            }

            // If we change from a situation of not contact toa  one of contact
            if (current_gap < 0.0 && previous_gap > 0.0) {
                double dynamic_factor = std::abs(current_gap)/std::abs(current_gap - previous_gap);
                dynamic_factor = (dynamic_factor > 1.0) ? 1.0 : dynamic_factor;
                KRATOS_DEBUG_ERROR_IF(dynamic_factor <= 0.0) << "DYNAMIC_FACTOR cannot be negative" << std::endl; // NOTE: THIS IS SUPPOSED TO BE IMPOSSIBLE (WE ARE USING ABS VALUES!!!!)
                rNode.SetValue(DYNAMIC_FACTOR, dynamic_factor);
            }
        }
    });

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

inline double ComputeDynamicFactorProcess::ComputeLogisticFactor(
    const double MaxGapThreshold,
    const double CurrentGap
    )
{
    const double exponent_factor = - 6.0 * (std::abs(CurrentGap)/MaxGapThreshold);
    return (1.0/(1.0 + std::exp(exponent_factor)));
}

} /// namespace Kratos
