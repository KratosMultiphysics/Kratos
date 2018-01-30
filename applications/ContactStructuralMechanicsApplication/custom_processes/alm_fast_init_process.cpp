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
        
        // We initialize the zero vector
        const array_1d<double, 3> zero_vector(3, 0.0);
        
        // We differentiate between frictional or frictionless
        const bool is_frictional = mrThisModelPart.Is(SLIP);
        
        // We initialize the penalty parameter
        const double& epsilon = mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY];
        
        bool init_delta_normal = false;
        Matrix zero_delta_normal;
        if (mrThisModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            init_delta_normal = true;
            const unsigned int dimension = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
            zero_delta_normal = ZeroMatrix( dimension, dimension );
        }
        
        // We iterate over the node
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for firstprivate(zero_vector)
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            // Weighted values
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
            if (is_frictional == true)
            {
                it_node->FastGetSolutionStepValue(WEIGHTED_SLIP) = 0.0;
            }
            
            // Penalty parameter
            it_node->SetValue(INITIAL_PENALTY, epsilon);
            
            // Nodal area
            it_node->SetValue(NODAL_AREA, 0.0);
            
            // Auxiliar values
            it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0);
            if (is_frictional == true)
            {
                it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, 0.0);
            }
            
            // The normal and tangents vectors
            it_node->SetValue(NORMAL, zero_vector);
            
            // The delta normal if necessary
            if (init_delta_normal == true)
            {
                it_node->SetValue(DELTA_NORMAL, zero_delta_normal);
            }
        }
        
        // Now we iterate over the conditions
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for firstprivate(zero_vector)
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            // The normal and tangents vectors
            it_cond->SetValue(NORMAL, zero_vector);
        }

        KRATOS_CATCH("");
    }
}
