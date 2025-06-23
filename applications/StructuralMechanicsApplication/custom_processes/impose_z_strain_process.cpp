// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// System includes
#include <algorithm>

// Project includes
#include "utilities/parallel_utilities.h"
#include "custom_processes/impose_z_strain_process.h"

namespace Kratos
{
ImposeZStrainProcess::ImposeZStrainProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeZStrainProcess::Execute()
{

}

/***********************************************************************************/
/***********************************************************************************/

void ImposeZStrainProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const ProcessInfo& r_proc_info = mrThisModelPart.GetProcessInfo();
    const double z_strain_value = mThisParameters["z_strain_value"].GetDouble();

    block_for_each(mrThisModelPart.Elements(), std::vector<double>(), [&r_proc_info, z_strain_value](
        auto& rElem, auto& rTls){
            const auto& r_int_pts = rElem.GetGeometry().IntegrationPoints(rElem.GetIntegrationMethod());
            if (r_int_pts.size() != rTls.size()){
                rTls.resize(r_int_pts.size());
                std::fill(rTls.begin(), rTls.end(), z_strain_value);
            }
            rElem.SetValuesOnIntegrationPoints(IMPOSED_Z_STRAIN_VALUE, rTls, r_proc_info);
    });

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ImposeZStrainProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name" : "please_specify_model_part_name",
        "z_strain_value"  : 0.01
    })" );

    return default_parameters;
}

// class ImposeZStrainProcess
} // namespace Kratos.
