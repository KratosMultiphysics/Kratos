 
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
#include "utilities/parallel_utilities.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/aalm_adapt_penalty_value_process.h"

namespace Kratos
{
void AALMAdaptPenaltyValueProcess::Execute()
{
    KRATOS_TRY;
    
    // We initialize the zero vector
    ProcessInfo& r_process_info = mrThisModelPart.GetProcessInfo();
    const double initial_penalty_parameter = r_process_info[INITIAL_PENALTY];
    const double max_gap_factor = r_process_info[MAX_GAP_FACTOR];
    const bool initialize = (r_process_info[STEP] == 1 && r_process_info[NL_ITERATION_NUMBER] == 1);
    
    // We iterate over the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    block_for_each(r_nodes_array, [&](Node& rNode) {
        // Initial value
        const double penalty_parameter = initialize ? initial_penalty_parameter : rNode.GetValue(INITIAL_PENALTY);
        
        // Weighted values
        const double nodal_area = rNode.Has(NODAL_AREA) ? rNode.GetValue(NODAL_AREA) : 0.0;
        if (nodal_area > std::numeric_limits<double>::epsilon()) {
            const double current_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
            const double previous_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP, 1)/nodal_area;
            
            // Nodal H x gap factor
            const double max_gap = max_gap_factor * rNode.FastGetSolutionStepValue(NODAL_H); // NOTE: This value must be studied
            
            // DEBUG
            KRATOS_TRACE("Current gap") << "current_gap: " << current_gap << "\tprevious_gap: " << previous_gap << "\tmax_gap: " << max_gap << std::endl;
            KRATOS_TRACE("Sign change gaps") << "(current_gap * previous_gap): " << (current_gap * previous_gap) << std::endl;
            KRATOS_TRACE("Compare previous gap") << "std::abs(previous_gap) > max_gap: " << (std::abs(previous_gap) > max_gap) << std::endl;
            KRATOS_TRACE("Compare current gap") << "std::abs(current_gap) > max_gap: " << (std::abs(current_gap) > max_gap) << std::endl;
            
            if ((current_gap * previous_gap) < 0.0) {
                if (std::abs(previous_gap) > max_gap) // NOTE: The abs is deduced from the paper (the algorithm is without abs)
                    rNode.SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (current_gap) * (std::abs(current_gap) + max_gap)/(current_gap - previous_gap)));
                else
                    rNode.SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (10.0 * current_gap)));
            } else if (std::abs(current_gap) > max_gap) { // NOTE: The abs is deduced from the paper (the algorithm is without abs)
                if (std::abs(current_gap - previous_gap) > std::max(current_gap/10.0, std::max(previous_gap/1.0, 5 * max_gap)))
                    rNode.SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter);
                else if ((std::abs(current_gap) <= std::abs(previous_gap) * 1.01 || std::abs(current_gap) >= std::abs(previous_gap) * 0.99) && (std::abs(current_gap) < 10.0 *  max_gap))
                    rNode.SetValue(INITIAL_PENALTY, penalty_parameter * std::pow((std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0), 2.0));
                else if (std::abs(current_gap) > std::abs(previous_gap) * 1.01)
                    rNode.SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter * (previous_gap/current_gap));
                else
                    rNode.SetValue(INITIAL_PENALTY, penalty_parameter * (std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0));
            } else
                rNode.SetValue(INITIAL_PENALTY, penalty_parameter);
        } else {
            rNode.SetValue(INITIAL_PENALTY, penalty_parameter);
        }
    });

    KRATOS_CATCH("");
} 
}
