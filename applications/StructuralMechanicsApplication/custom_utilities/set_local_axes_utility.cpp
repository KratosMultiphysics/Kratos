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
// Project includes
#include "set_local_axes_utility.h"

namespace Kratos {

void SetLocalAxesUtility::SetLocalAxisCartesianSystem(
          ModelPart &rModelPart,
          Parameters ThisParameters)
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

} // namespace Kratos