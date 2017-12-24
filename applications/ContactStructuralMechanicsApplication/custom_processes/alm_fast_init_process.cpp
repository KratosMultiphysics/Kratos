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
#include "custom_processes/alm_fast_init_process.h"

namespace Kratos
{
void ALMFastInit::Execute()
{
    KRATOS_TRY;
    
    // We differentiate between frictional or frictionless
    const bool is_frictional = mrThisModelPart.Is(SLIP);
    
    // We initialize the penalty parameter
    const double epsilon = mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY];
    
    // We iterate over the nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
        // Weighted values
        it_node->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
        if (is_frictional == true)
            it_node->FastGetSolutionStepValue(WEIGHTED_SLIP) = 0.0;
        
        // Penalty parameter
        it_node->SetValue(INITIAL_PENALTY, epsilon);
        
        // Auxiliar values
        it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0);
        if (is_frictional == true)
            it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, 0.0);
    }
    
    // Now we iterate over the conditions
    ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)
        (conditions_array.begin() + i)->SetValue(NORMAL, ZeroVector(3)); // The normal and tangents vectors


    KRATOS_CATCH("");
} // class ALMFastInit
} // namespace Kratos
