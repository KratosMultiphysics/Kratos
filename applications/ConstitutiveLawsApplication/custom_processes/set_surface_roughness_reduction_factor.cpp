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
#include "custom_processes/set_surface_roughness_reduction_factor.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
SetSurfaceRoughnessReductionFactor::SetSurfaceRoughnessReductionFactor(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetSurfaceRoughnessReductionFactor::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY
    const double surface_roughness = mThisParameters["surface_roughness"].GetDouble();
    const double mat_param_ar = mThisParameters["mat_param_ar"].GetDouble();
    const double mat_param_sigmautsmin = mThisParameters["mat_param_sigmautsmin"].GetDouble();

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        rElement.SetValue(SURFACE_ROUGHNESS, surface_roughness);
        rElement.SetValue(MATERIAL_PARAMETER_C1, mat_param_ar);
        rElement.SetValue(MATERIAL_PARAMETER_C2, mat_param_sigmautsmin);
    });
        
        KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetSurfaceRoughnessReductionFactor::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "surface_roughness"         : 1.0,
        "mat_param_ar"              : 0.0,
        "mat_param_sigmautsmin"     : 1.0
    })");

    return default_parameters;
}

} // namespace Kratos.