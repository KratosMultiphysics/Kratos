//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <tuple>

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Include base h
#include "rans_adjoint_utilities.h"

namespace Kratos
{
double AdjointUtilities::CalculateWallHeightDerivative(
    const ConditionType& rCondition,
    const array_1d<double, 3>& rNormal,
    const array_1d<double, 3>& rNormalDerivative)
{
    KRATOS_TRY

    // for some weird reason, I cannot make the following array_1d<double, 3> to auto.
    // Clang is compiling fine, and it works as it suppose to be. But in gcc, it compiles
    // but all the tests start to fail not by crashing, but giving false values.
    const array_1d<double, 3>& normal = rNormal / norm_2(rNormal);
    const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];

    const auto& r_parent_geometry = r_parent_element.GetGeometry();
    const auto& r_condition_geometry = rCondition.GetGeometry();

    const auto& parent_center = r_parent_geometry.Center();
    const auto& condition_center = r_condition_geometry.Center();

    return inner_prod(condition_center - parent_center, normal);

    KRATOS_CATCH("");
}

void AdjointUtilities::CalculateYPlusAndUtauDerivative(
    double& rYPlusDerivative,
    double& rUTauDerivative,
    const double YPlus,
    const double UTau,
    const double WallVelocity,
    const double WallVelocityDerivative,
    const double WallHeight,
    const double WallHeightDerivative,
    const double KinematicViscosity,
    const double Kappa,
    const double Beta,
    const double YPlusLimit)
{
    KRATOS_TRY

    if (YPlus > YPlusLimit) {
        // compute logarithmic wall law derivatives
        rUTauDerivative = (WallVelocityDerivative -
                           UTau * WallHeightDerivative / (Kappa * WallHeight)) /
                          (1 / Kappa + WallVelocity / UTau);
    } else {
        // compute linear wall law derivatives
        rUTauDerivative =
            (WallVelocityDerivative - std::pow(UTau, 2) * WallHeightDerivative / KinematicViscosity) /
            (UTau * (WallHeight / KinematicViscosity + WallVelocity / std::pow(UTau, 2)));
    }

    rYPlusDerivative =
        (rUTauDerivative * WallHeight + UTau * WallHeightDerivative) / KinematicViscosity;

    KRATOS_CATCH("");
}

} // namespace Kratos