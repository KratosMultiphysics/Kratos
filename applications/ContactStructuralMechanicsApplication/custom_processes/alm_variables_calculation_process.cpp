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
#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_utilities/contact_utilities.h"
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
    void ALMVariablesCalculationProcess::Execute()
    {
        KRATOS_TRY
        
        /* We compute the penalty factor */
        
        // We initialize the mean values
        double mean_young_modulus_slave  = 0.0;
        double mean_nodal_h_slave        = 0.0;
        double mean_young_modulus_master = 0.0;
        double mean_nodal_h_master       = 0.0;
        
        // We initialize the total areas and volumes
        double total_volume_slave  = 0.0;
        double total_area_slave    = 0.0;
        double total_volume_master = 0.0;
        double total_area_master   = 0.0;
        
        // Now we iterate over the conditions to calculate the nodal area
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for reduction(+:total_volume_slave, total_area_slave, mean_young_modulus_slave, mean_nodal_h_slave, total_volume_master, total_area_master, mean_young_modulus_master, mean_nodal_h_master)
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            // We get the condition geometry
            GeometryType& r_this_geometry = it_cond->GetGeometry();
            const unsigned num_nodes_geometry = r_this_geometry.size();
            
            // We get the values from the condition
            Kratos::Properties& this_properties = it_cond->GetProperties();
            const double& young_modulus = this_properties[YOUNG_MODULUS];
            const double& element_volume = it_cond->GetGeometry().Area();
            
            // We get the values from the condition
            const double& condition_area = r_this_geometry.Area();
            const double nodal_condition_area = condition_area/num_nodes_geometry;
            
            if (it_cond->Is(SLAVE) == true)
            {
                total_volume_slave += element_volume;
                total_area_slave += condition_area;
                mean_young_modulus_slave += young_modulus * element_volume;
                
                for (unsigned int i_node = 0; i_node < num_nodes_geometry; i_node++)
                {
                    mean_nodal_h_slave += r_this_geometry[i_node].FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                }
            }
            
            if (it_cond->Is(MASTER) == true)
            {
                total_volume_master += element_volume;
                total_area_master += condition_area;
                mean_young_modulus_master += young_modulus * element_volume;
                
                for (unsigned int i_node = 0; i_node < num_nodes_geometry; i_node++)
                {
                    mean_nodal_h_master += r_this_geometry[i_node].FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                }
            }
        }
        
        // Now we divide between the total areas and volumes
        mean_nodal_h_slave /= (total_area_slave + 1.0e-12); 
        mean_young_modulus_slave /= (total_volume_slave + 1.0e-12);
        
        mean_nodal_h_master /= (total_area_master  + 1.0e-12); 
        mean_young_modulus_master /= (total_volume_master + 1.0e-12);
        
        // Finally we compute the penalty factor
        const double penalty_parameter_slave  = mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + 1.0e-12);
        const double scale_factor_slave    = mPenaltyScale * mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + 1.0e-12);
        const double penalty_parameter_master = mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + 1.0e-12);
        const double scale_factor_master   = mPenaltyScale * mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + 1.0e-12); 
        
        mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY] = (penalty_parameter_slave > penalty_parameter_master) ? penalty_parameter_slave : penalty_parameter_master; // NOTE: > or <? , we are supposed to take the largest of the values (more stiff)
        mrThisModelPart.GetProcessInfo()[SCALE_FACTOR] = (scale_factor_slave > scale_factor_master) ? scale_factor_slave : scale_factor_master;
        
        KRATOS_CATCH("")
    }
}
