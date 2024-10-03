// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
//
//  Main authors:    Luis Antonio Goncalves Junior
//                   Alejandro Cornejo
//

// System includes

// External includes



// Project includes
#include "includes/model_part.h"
#include "custom_processes/set_punching_stress_concentration_factor.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
SetPunchingStressConcentrationFactor::SetPunchingStressConcentrationFactor(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetPunchingStressConcentrationFactor::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY
    const double hole_diameter = mThisParameters["hole_diameter"].GetDouble();
    const double specimen_width = mThisParameters["specimen_width"].GetDouble();

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        rElement.SetValue(HOLE_DIAMETER, hole_diameter);
        rElement.SetValue(SPECIMEN_WIDTH, specimen_width);
    });
        
        KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetPunchingStressConcentrationFactor::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "hole_diameter"            : 0.0,
        "specimen_width"           : 1.0
    })");

    return default_parameters;
}

} // namespace Kratos.