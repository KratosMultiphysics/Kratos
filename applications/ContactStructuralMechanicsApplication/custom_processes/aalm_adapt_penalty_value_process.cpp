 
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

// System includes

// External includes

// Project includes
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
    
    // We iterate over the node
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        
        // Initial value
        const double penalty_parameter = initialize ? initial_penalty_parameter : it_node->GetValue(INITIAL_PENALTY);
        
        // Weighted values
        const double nodal_area = it_node->Has(NODAL_AREA) ? it_node->GetValue(NODAL_AREA) : 0.0;
        if (nodal_area > std::numeric_limits<double>::epsilon()) {
            const double current_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
            const double previous_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1)/nodal_area;
            
            // Nodal H x gap factor
            const double max_gap = max_gap_factor * it_node->FastGetSolutionStepValue(NODAL_H); // NOTE: This value must be studied
            
            // DEBUG
            KRATOS_TRACE("Current gap") << "current_gap: " << current_gap << "\tprevious_gap: " << previous_gap << "\tmax_gap: " << max_gap << std::endl;
            KRATOS_TRACE("Sign change gaps") << "(current_gap * previous_gap): " << (current_gap * previous_gap) << std::endl;
            KRATOS_TRACE("Compare previous gap") << "std::abs(previous_gap) > max_gap: " << (std::abs(previous_gap) > max_gap) << std::endl;
            KRATOS_TRACE("Compare current gap") << "std::abs(current_gap) > max_gap: " << (std::abs(current_gap) > max_gap) << std::endl;
            
            if ((current_gap * previous_gap) < 0.0) {
                if (std::abs(previous_gap) > max_gap) // NOTE: The abs is deduced from the paper (the algorithm is without abs)
                    it_node->SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (current_gap) * (std::abs(current_gap) + max_gap)/(current_gap - previous_gap)));
                else
                    it_node->SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (10.0 * current_gap)));
            } else if (std::abs(current_gap) > max_gap) { // NOTE: The abs is deduced from the paper (the algorithm is without abs)
                if (std::abs(current_gap - previous_gap) > std::max(current_gap/10.0, std::max(previous_gap/1.0, 5 * max_gap)))
                    it_node->SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter);
                else if ((std::abs(current_gap) <= std::abs(previous_gap) * 1.01 || std::abs(current_gap) >= std::abs(previous_gap) * 0.99) && (std::abs(current_gap) < 10.0 *  max_gap))
                    it_node->SetValue(INITIAL_PENALTY, penalty_parameter * std::pow((std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0), 2.0));
                else if (std::abs(current_gap) > std::abs(previous_gap) * 1.01)
                    it_node->SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter * (previous_gap/current_gap));
                else
                    it_node->SetValue(INITIAL_PENALTY, penalty_parameter * (std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0));
            } else
                it_node->SetValue(INITIAL_PENALTY, penalty_parameter);
            } else {
                it_node->SetValue(INITIAL_PENALTY, penalty_parameter);
        }
    }

    KRATOS_CATCH("");
} 
}
