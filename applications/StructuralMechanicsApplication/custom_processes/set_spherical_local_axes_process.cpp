// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#include "includes/model_part.h"
#include "custom_processes/set_spherical_local_axes_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetSphericalLocalAxesProcess::SetSphericalLocalAxesProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetSphericalLocalAxesProcess::ExecuteInitialize()
{
    KRATOS_TRY
    const array_1d<double, 3>& spherical_reference_axis  = mThisParameters["spherical_reference_axis"].GetVector();
    const array_1d<double, 3>& spherical_central_point   = mThisParameters["spherical_central_point"].GetVector();
    const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(MathUtils<double>::Norm3(spherical_reference_axis) < tolerance) << "The spherical_reference_axis has norm zero" << std::endl;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        array_1d<double, 3> local_axis_1;

        const array_1d<double, 3> coords = rElement.GetGeometry().Center();
        noalias(local_axis_1) = coords - spherical_central_point;
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
        rElement.SetValue(LOCAL_AXIS_1, local_axis_1);

        if (mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
            array_1d<double, 3> local_axis_2;
            if (std::abs(local_axis_1[2]) <= tolerance ||
                std::abs(local_axis_1[1] * spherical_reference_axis[2] - local_axis_1[2] * spherical_reference_axis[1]) <= tolerance) {
                noalias(local_axis_2) = spherical_reference_axis;
                ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);
            }
            // Let's compute the second local axis
            const double x = 1.0;
            const double y = -(local_axis_1[0] * spherical_reference_axis[2] - local_axis_1[2] * spherical_reference_axis[0]) / (local_axis_1[1] * spherical_reference_axis[2] - local_axis_1[2] * spherical_reference_axis[1]);
            const double z = -local_axis_1[0] / local_axis_1[2] - local_axis_1[1] / local_axis_1[2] * y;
            local_axis_2[0] = x;
            local_axis_2[1] = y;
            local_axis_2[2] = z;
            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);
            rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        }
    });
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetSphericalLocalAxesProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "spherical_reference_axis"   : [0.0,0.0,1.0],
        "spherical_central_point"    : [0.0,0.0,0.0]
    })" );

    return default_parameters;
}

// class SetSphericalLocalAxesProcess
} // namespace Kratos.
