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

// System includes

// External includes
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
// Project includes
#include "set_local_axes_utility.h"

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
    BoundedVector<double, 3> local_axis_1;
    BoundedVector<double, 3> local_axis_2;
    BoundedVector<double, 3> local_axis_3;

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

    const BoundedVector<double, 3> generatrix_axis  = ThisParameters["cylindrical_generatrix_axis"].GetVector();
    const BoundedVector<double, 3> generatrix_point = ThisParameters["cylindrical_generatrix_point"].GetVector();

    BoundedVector<double, 3> local_axis_1;
    BoundedVector<double, 3> local_axis_2;
    BoundedVector<double, 3> local_axis_3;

    block_for_each(rModelPart.Elements(), [&](Element &rElement) {
        const BoundedVector<double, 3> coords = rElement.GetGeometry().Center();
        const double c = -generatrix_axis(0) * coords(0) - generatrix_axis(1) * coords(1) - generatrix_axis(2) * coords(2);
        const double lambda = -(generatrix_axis(0) * generatrix_point(0) + generatrix_axis(1) * generatrix_point(1) + generatrix_axis(2) * generatrix_point(2) + c) / (std::pow(generatrix_axis(0), 2) + std::pow(generatrix_axis(1), 2) + std::pow(generatrix_axis(2), 2));

        BoundedVector<double, 3> intersection;
        noalias(intersection) = generatrix_point + lambda * generatrix_axis;

        noalias(local_axis_1) = coords - intersection;
        noalias(local_axis_2) = generatrix_axis;
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

} // namespace Kratos