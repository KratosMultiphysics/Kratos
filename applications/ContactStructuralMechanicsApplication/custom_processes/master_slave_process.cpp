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
#include "custom_processes/master_slave_process.h"

namespace Kratos
{
    void MasterSlaveProcess::Execute()
    {
        KRATOS_TRY;
        
        // Now we iterate over the conditions
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            const auto& this_geometry = it_cond->GetGeometry();
            
            // We set the condition as master or slave (master by default)
            bool is_slave = true;
            for (unsigned int i_node = 0; i_node < this_geometry.size(); ++i_node)
            {
                if (!this_geometry[i_node].Is(SLAVE)) is_slave = false;
            }
            if (is_slave == true)  it_cond->Set(SLAVE, true);
            else it_cond->Set(MASTER, true);
        }

        KRATOS_CATCH("");
    }
}
