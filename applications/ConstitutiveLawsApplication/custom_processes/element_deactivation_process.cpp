// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes



// Project includes
#include "includes/model_part.h"
#include "custom_processes/element_deactivation_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
ElementDeactivationProcess::ElementDeactivationProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ) : mrThisModelPart(rThisModelPart),
    mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void ElementDeactivationProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY








    KRATOS_CATCH("")
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ElementDeactivationProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "variable_name"               : "DAMAGE",
        "variable_maximum_threshold"  : 0.9999,
        "average_calculation_over_ip" : true
    })");

    return default_parameters;
}

} // namespace Kratos.