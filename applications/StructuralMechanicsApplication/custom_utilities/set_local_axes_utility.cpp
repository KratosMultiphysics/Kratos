// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"

// Project includes
#include "custom_utilities/set_local_axes_utility.h"

namespace Kratos {

void SetLocalAxesUtility::SetLocalAxisCartesianSystem(
    ModelPart &rModelPart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "local_axes_coordinate_system"  : "cartesian",
        "cartesian_local_axis"          : [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    })");
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    const Matrix cartesian_local_axes_matrix = ThisParameters["cartesian_local_axis"].GetMatrix();
    array_1d<double, 3> local_axis_1;
    array_1d<double, 3> local_axis_2;
    array_1d<double, 3> local_axis_3;

    local_axis_1(0) = cartesian_local_axes_matrix(0, 0);
    local_axis_1(1) = cartesian_local_axes_matrix(0, 1);
    local_axis_1(2) = cartesian_local_axes_matrix(0, 2);

    local_axis_2(0) = cartesian_local_axes_matrix(1, 0);
    local_axis_2(1) = cartesian_local_axes_matrix(1, 1);
    local_axis_2(2) = cartesian_local_axes_matrix(1, 2);

    local_axis_3(0) = cartesian_local_axes_matrix(2, 0);
    local_axis_3(1) = cartesian_local_axes_matrix(2, 1);
    local_axis_3(2) = cartesian_local_axes_matrix(2, 2);

    CheckAndNormalizeVector(local_axis_1);
    CheckAndNormalizeVector(local_axis_2);
    CheckAndNormalizeVector(local_axis_3);

    block_for_each(rModelPart.Elements(), [&](Element& rElement) {
        rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
        rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        rElement.SetValue(LOCAL_AXIS_3, local_axis_3);
    });

    KRATOS_CATCH("")
}

void SetLocalAxesUtility::SetLocalAxisCylindricalSystem(
    ModelPart &rModelPart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "local_axes_coordinate_system"  : "cylindrical",
        "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
        "cylindrical_generatrix_point"  : [0.0,0.0,0.0]
    })");
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    const array_1d<double, 3>& r_generatrix_axis  = ThisParameters["cylindrical_generatrix_axis"].GetVector();
    const array_1d<double, 3>& r_generatrix_point = ThisParameters["cylindrical_generatrix_point"].GetVector();

    KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;

    block_for_each(rModelPart.Elements(), [&](Element &rElement) {
        array_1d<double, 3> local_axis_1;
        array_1d<double, 3> local_axis_2;
        array_1d<double, 3> local_axis_3;
        const array_1d<double, 3> coords = rElement.GetGeometry().Center();
        const double c = -r_generatrix_axis(0) * coords(0) - r_generatrix_axis(1) * coords(1) - r_generatrix_axis(2) * coords(2);
        const double lambda = -(r_generatrix_axis(0) * r_generatrix_point(0) + r_generatrix_axis(1) * r_generatrix_point(1) + r_generatrix_axis(2) * r_generatrix_point(2) + c) / (std::pow(r_generatrix_axis(0), 2) + std::pow(r_generatrix_axis(1), 2) + std::pow(r_generatrix_axis(2), 2));

        array_1d<double, 3> intersection;
        noalias(intersection) = r_generatrix_point + lambda * r_generatrix_axis;

        noalias(local_axis_1) = coords - intersection;
        noalias(local_axis_2) = r_generatrix_axis;
        noalias(local_axis_3) = MathUtils<double>::CrossProduct(local_axis_1, local_axis_2);

        CheckAndNormalizeVector(local_axis_1);
        CheckAndNormalizeVector(local_axis_2);
        CheckAndNormalizeVector(local_axis_3);

        rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
        rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        rElement.SetValue(LOCAL_AXIS_3, local_axis_3);
    });

    KRATOS_CATCH("")
}

void SetLocalAxesUtility::SetLocalAxisSphericalSystem(
    ModelPart &rModelPart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "local_axes_coordinate_system"  : "spherical",
        "spherical_reference_axis"   : [0.0,0.0,1.0],
        "spherical_central_point"    : [0.0,0.0,0.0]
    })");
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    const array_1d<double, 3> spherical_reference_axis  = ThisParameters["spherical_reference_axis"].GetVector();
    const array_1d<double, 3> spherical_central_point   = ThisParameters["spherical_central_point"].GetVector();
    const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(MathUtils<double>::Norm3(spherical_reference_axis) < tolerance) << "The spherical_reference_axis has norm zero" << std::endl;

    block_for_each(rModelPart.Elements(), [&](Element &rElement) {
        array_1d<double, 3> local_axis_1;
        array_1d<double, 3> local_axis_2;
        array_1d<double, 3> local_axis_3;

        const array_1d<double, 3> coords = rElement.GetGeometry().Center();
        noalias(local_axis_1) = coords - spherical_central_point;
        CheckAndNormalizeVector(local_axis_1);

        if (std::abs(local_axis_1(2)) <= tolerance ||
            std::abs(local_axis_1(1) * spherical_reference_axis(2) - local_axis_1(2) * spherical_reference_axis(1)) <= tolerance) {
            noalias(local_axis_2) = spherical_reference_axis;
            CheckAndNormalizeVector(local_axis_2);
            noalias(local_axis_3) = MathUtils<double>::CrossProduct(local_axis_1, local_axis_2);
            CheckAndNormalizeVector(local_axis_3);
        }
        // Let's compute the second local axis
        const double x = 1.0;
        const double y = -(local_axis_1(0) * spherical_reference_axis(2) - local_axis_1(2) * spherical_reference_axis(0)) / (local_axis_1(1) * spherical_reference_axis(2) - local_axis_1(2) * spherical_reference_axis(1));
        const double z = -local_axis_1(0) / local_axis_1(2) - local_axis_1(1) / local_axis_1(2) * y;
        local_axis_2(0) = x;
        local_axis_2(1) = y;
        local_axis_2(2) = z;
        CheckAndNormalizeVector(local_axis_2);

        noalias(local_axis_3) = MathUtils<double>::CrossProduct(local_axis_1, local_axis_2);
        CheckAndNormalizeVector(local_axis_3);

        rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
        rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        rElement.SetValue(LOCAL_AXIS_3, local_axis_3);
    });

    KRATOS_CATCH("")
}

} // namespace Kratos
