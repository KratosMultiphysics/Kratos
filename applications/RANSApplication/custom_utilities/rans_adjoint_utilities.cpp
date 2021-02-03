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

double RansAdjointUtilities::CalculateVectorNormDerivative(
    const double VectorNorm,
    const array_1d<double, 3>& rVector,
    const array_1d<double, 3>& rVectorDerivative)
{
    if (VectorNorm > 0.0) {
        return inner_prod(rVector, rVectorDerivative) / VectorNorm;
    } else {
        return 0.0;
    }
}

array_1d<double, 3> RansAdjointUtilities::CalculateUnitVectorDerivative(
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

double RansAdjointUtilities::CalculateWallHeightConditionDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentElementGeometry.Center();

    array_1d<double, 3> condition_center_derivative = ZeroVector(3);
    condition_center_derivative[DirectionIndex] = 1.0 / rConditionGeometry.PointsNumber();

    array_1d<double, 3> parent_center_derivative = ZeroVector(3);
    parent_center_derivative[DirectionIndex] = 1.0 / rParentElementGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) +
           inner_prod(condition_center_derivative - parent_center_derivative, rUnitNormal);
}

double RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentElementGeometry.Center();

    array_1d<double, 3> parent_center_derivative = ZeroVector(3);
    parent_center_derivative[DirectionIndex] = 1.0 / rParentElementGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) -
           inner_prod(parent_center_derivative, rUnitNormal);
}

void RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
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

template<>
double RansAdjointUtilities::GeometricalDerivatives<2, 2>::DomainSizeDerivative(
    const GeometryType& rGeometry,
    const IndexType NodeIndex,
    const IndexType DirectionIndex)
{
    const double lx = rGeometry[0].X() - rGeometry[1].X();
    const double lx_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 0);

    const double ly = rGeometry[0].Y() - rGeometry[1].Y();
    const double ly_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 1);

    const double length = lx * lx + ly * ly;
    const double length_derivative = 2.0 * lx * lx_derivative + 2.0 * ly * ly_derivative;

    const double domain_size = std::sqrt(length);
    const double domain_size_derivative = 0.5 * length_derivative / std::sqrt(length);

    return domain_size_derivative;
}

template<>
double RansAdjointUtilities::GeometricalDerivatives<3, 3>::DomainSizeDerivative(
    const GeometryType& rGeometry,
    const IndexType NodeIndex,
    const IndexType DirectionIndex)
{
    // const double a = MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1));
    // // const double a_derivative =

    // const double b = MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2));
    // const double c = MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0));

    // const double s = (a+b+c) / 2.0;

    // return std::sqrt(s*(s-a)*(s-b)*(s-c));
}

} // namespace Kratos