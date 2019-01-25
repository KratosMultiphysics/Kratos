
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

    // Getting delta time
    const double delta_time = r_process_info[DELTA_TIME];

    // Getting logistic factor
    double logistic_factor = 1.0;

    // Impact time duration
    const double impact_time_duration = r_process_info.Has(IMPACT_TIME_DURATION) ? r_process_info.GetValue(IMPACT_TIME_DURATION) : 0.0;

    // We iterate over the node
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    #pragma omp parallel for firstprivate(logistic_factor)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        // Computing only on SLAVE nodes
        if (it_node->Is(SLAVE) && it_node->Is(ACTIVE)) {
            // Weighted values
            const double current_gap  = it_node->FastGetSolutionStepValue(WEIGHTED_GAP);
            const double previous_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1);

            // Reseting the logistic factor
            logistic_factor = 1.0;

            // Computing actual logistic factor
            if (impact_time_duration > 0.0 && current_gap <= 0.0) {
                double& node_delta_time = it_node->GetValue(DELTA_TIME);
                node_delta_time += delta_time;
                logistic_factor = ComputeLogisticFactor(impact_time_duration, node_delta_time);
            } else if (impact_time_duration > 0.0 && current_gap > 0.0) {
                it_node->SetValue(DELTA_TIME, 0.0);
            }

            // If we change from a situation of not contact toa  one of contact
            if (current_gap < 0.0 && previous_gap > 0.0) {
                double dynamic_factor = std::abs(current_gap)/std::abs(current_gap - previous_gap);
                dynamic_factor = (dynamic_factor > 1.0) ? 1.0 : dynamic_factor;
                KRATOS_DEBUG_ERROR_IF(dynamic_factor <= 0.0) << "DYNAMIC_FACTOR cannot be negative" << std::endl; // NOTE: THIS IS SUPPOSED TO BE IMPOSSIBLE (WE ARE USING ABS VALUES!!!!)
                it_node->SetValue(DYNAMIC_FACTOR, logistic_factor * dynamic_factor);
            } else {
                it_node->SetValue(DYNAMIC_FACTOR, logistic_factor);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeDynamicFactorProcess::ExecuteInitialize()
{
    // Getting process info
    ProcessInfo& r_process_info = mrThisModelPart.GetProcessInfo();

    // We initialize in case IMPACT_TIME_DURATION
    if (r_process_info.Has(IMPACT_TIME_DURATION)) {
        // We initialize the DELTA_TIME on the nodes
        VariableUtils().SetNonHistoricalVariable(DELTA_TIME, 0.0, mrThisModelPart.Nodes());
    }
}

/***********************************************************************************/
/***********************************************************************************/

double ComputeDynamicFactorProcess::ComputeLogisticFactor(
    const double ImpactTimeDuration,
    const double CurrentDeltaTime
    )
{
    const double exponent_factor = - 6.0 * (CurrentDeltaTime/ImpactTimeDuration);
    return (1.0/(1.0 + std::exp(exponent_factor)));
}

} /// namespace Kratos
