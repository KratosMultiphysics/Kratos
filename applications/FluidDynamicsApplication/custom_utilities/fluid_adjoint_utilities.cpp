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

// External includes

// Project includes
#include "includes/define.h"

// Application includes

// Include base h
#include "fluid_adjoint_utilities.h"

namespace Kratos
{

/***************************************************************************************/
/*************************************** double ****************************************/
/***************************************************************************************/

template <>
template <>
std::array<const Variable<double>*, 2> FluidAdjointUtilities<2>::GetRelevantGradientVariableComponentList<double, 3>(
    const IndexType DirectionIndex,
    const Variable<double>& rVariable,
    const std::array<const Variable<double>*, 3>& rAllGradientVariableComponents)
{
    return {rAllGradientVariableComponents[0], rAllGradientVariableComponents[1]};
}

template <>
template <>
std::array<const Variable<double>*, 3> FluidAdjointUtilities<3>::GetRelevantGradientVariableComponentList<double, 3>(
    const IndexType DirectionIndex,
    const Variable<double>& rVariable,
    const std::array<const Variable<double>*, 3>& rAllGradientVariableComponents)
{
    return {rAllGradientVariableComponents[0], rAllGradientVariableComponents[1], rAllGradientVariableComponents[2]};
}

/***************************************************************************************/
/********************************* array_1d<double, 3> *********************************/
/***************************************************************************************/

template <>
template <>
const Variable<double>& FluidAdjointUtilities<2>::GetRelevantVariable<array_1d<double, 3>>(
    const IndexType DirectionIndex,
    const Variable<array_1d<double, 3>>& rVariable,
    const std::array<const Variable<double>*, 3>& rAllVariableComponents)
{
    return *rAllVariableComponents[DirectionIndex];
}

template <>
template <>
const Variable<double>& FluidAdjointUtilities<3>::GetRelevantVariable<array_1d<double, 3>>(
    const IndexType DirectionIndex,
    const Variable<array_1d<double, 3>>& rVariable,
    const std::array<const Variable<double>*, 3>& rAllVariableComponents)
{
    return *rAllVariableComponents[DirectionIndex];
}

template <>
template <>
std::array<const Variable<double>*, 2> FluidAdjointUtilities<2>::GetRelevantGradientVariableComponentList<array_1d<double, 3>, 9>(
    const IndexType DirectionIndex,
    const Variable<array_1d<double, 3>>& rVariable,
    const std::array<const Variable<double>*, 9>& rAllGradientVariableComponents)
{
    return {
        rAllGradientVariableComponents[DirectionIndex * 3],
        rAllGradientVariableComponents[DirectionIndex * 3 + 1]};
}

template<>
template<>
std::array<const Variable<double>*, 3> FluidAdjointUtilities<3>::GetRelevantGradientVariableComponentList<array_1d<double, 3>, 9>(
    const IndexType DirectionIndex,
    const Variable<array_1d<double, 3>>& rVariable,
    const std::array<const Variable<double>*, 9>& rAllGradientVariableComponents)
{
    return {
        rAllGradientVariableComponents[DirectionIndex * 3],
        rAllGradientVariableComponents[DirectionIndex * 3 + 1],
        rAllGradientVariableComponents[DirectionIndex * 3 + 2]};
}

template<>
double FluidAdjointUtilities<2>::CalculateTriangleAreaDerivative(
    const GeometryType& rGeometry,
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex)
{
    KRATOS_TRY

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    const double detJ_derivative = x10 * y20_derivative + x10_derivative * y20 - y10 * x20_derivative - y10_derivative * x20;

    return 0.5 * detJ_derivative;

    KRATOS_CATCH("");
}

template<>
double FluidAdjointUtilities<3>::CalculateTriangleAreaDerivative(
    const GeometryType& rGeometry,
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex)
{
    KRATOS_TRY

    const array_1d<double, 3>& l01 = rGeometry[0].Coordinates() - rGeometry[1].Coordinates();
    const array_1d<double, 3> l01_derivative{
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 0, 1, 0),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 0, 1, 1),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 0, 1, 2)
    };
    const double a = norm_2(l01);
    const double a_derivative = inner_prod(l01, l01_derivative) / a;

    const array_1d<double, 3>& l12 = rGeometry[1].Coordinates() - rGeometry[2].Coordinates();
    const array_1d<double, 3> l12_derivative{
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 2, 0),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 2, 1),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 2, 2)
    };
    const double b = norm_2(l12);
    const double b_derivative = inner_prod(l12, l12_derivative) / b;

    const array_1d<double, 3>& l20 = rGeometry[2].Coordinates() - rGeometry[0].Coordinates();
    const array_1d<double, 3> l20_derivative{
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1),
        EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 2)
    };
    const double c = norm_2(l20);
    const double c_derivative = inner_prod(l20, l20_derivative) / c;

    const double s = (a + b + c) / 2.0;
    const double s_derivative = (a_derivative + b_derivative + c_derivative) / 2.0;

    const double area = std::sqrt(s * (s - a) * (s - b) * (s - c));

    double area_derivative = 0;
    area_derivative += s_derivative * (s - a) * (s - b) * (s - c);
    area_derivative += s * (s_derivative - a_derivative) * (s - b) * (s - c);
    area_derivative += s * (s - a) * (s_derivative - b_derivative) * (s - c);
    area_derivative += s * (s - a) * (s - b) * (s_derivative - c_derivative);
    area_derivative *= 0.5 / area;

    return area_derivative;

    KRATOS_CATCH("");
}

// template instantiations

template class FluidAdjointUtilities<2>;
template class FluidAdjointUtilities<3>;

} // namespace Kratos