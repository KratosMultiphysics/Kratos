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

// template instantiations

template class FluidAdjointUtilities<2>;
template class FluidAdjointUtilities<3>;

} // namespace Kratos