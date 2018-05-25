
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
#include "custom_processes/compute_dynamic_factor_process.h"

namespace Kratos
{
void ComputeDynamicFactorProcess::Execute()
{
    KRATOS_TRY;

    // We initialize the zero vector
//     ProcessInfo& r_process_info = mrThisModelPart.GetProcessInfo();
    const double distance_threshold = 0.0;
//     const double distance_threshold = r_process_info[DISTANCE_THRESHOLD]; // TODO: Think about this!!!!!

    // We iterate over the node
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;

        // Weighted values
        const double nodal_area = it_node->GetValue(NODAL_AREA);
        const double current_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
        const double previous_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1)/nodal_area;

        // If we change from a situation of not contact toa  one of contact
        if (current_gap < distance_threshold && previous_gap > previous_gap) {
            double dynamic_factor = std::abs(current_gap - distance_threshold)/std::abs(current_gap - previous_gap);
            dynamic_factor = (dynamic_factor > 1.0) ? 1.0 :dynamic_factor;
            KRATOS_DEBUG_ERROR_IF(dynamic_factor <= 0.0) << "DYNAMIC_FACTOR cannot be negative" << std::endl; // NOTE: THIS IS SUPPOSED TO BE IMPOSSIBLE (WE ARE USING ABS VALUES!!!!)
            it_node->SetValue(DYNAMIC_FACTOR, dynamic_factor);
        } else
            it_node->SetValue(DYNAMIC_FACTOR, 1.0);
    }

    KRATOS_CATCH("");
}
}
