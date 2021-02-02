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

array_1d<double, 3> AdjointUtilities::CalculateUnitVectorDerivative(
    const double VectorMagnitude,
    const array_1d<double, 3>& rUnitVector,
    const array_1d<double, 3>& rVectorDerivative)
{
    if (VectorMagnitude > 0.0) {
        return (rVectorDerivative - rUnitVector * inner_prod(rUnitVector, rVectorDerivative)) /
               VectorMagnitude;
    } else {
        return ZeroVector(3);
    }
}

double AdjointUtilities::CalculateWallHeightConditionDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentGeometry.Center();

    array_1d<double, 3> condition_center_derivative = ZeroVector(3);
    condition_center_derivative[DirectionIndex] = 1.0 / rConditionGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) +
           inner_prod(condition_center_derivative, rUnitNormal);
}

double AdjointUtilities::CalculateWallHeightParentElementDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentGeometry.Center();

    array_1d<double, 3> parent_center_derivative = ZeroVector(3);
    parent_center_derivative[DirectionIndex] = 1.0 / rParentGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) -
           inner_prod(parent_center_derivative, rUnitNormal);
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