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

// External includes

// Project includes

// Include base h
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{

template<>
void FluidCalculationUtilities::AssignValue(
    array_1d<double, 2>& rOutput,
    const array_1d<double, 3>& rInput)
{
    rOutput[0] = rInput[0];
    rOutput[1] = rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::AssignValue(
    TOutputDataType& rOutput,
    const TInputDataType& rInput)
{
    rOutput = rInput;
}

template<>
void FluidCalculationUtilities::UpdateValue(
    double& rOutput,
    const double& rInput)
{
    rOutput += rInput;
}

template<>
void FluidCalculationUtilities::UpdateValue(
    array_1d<double, 2>& rOutput,
    const array_1d<double, 3>& rInput)
{
    rOutput[0] += rInput[0];
    rOutput[1] += rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::UpdateValue(
    TOutputDataType& rOutput,
    const TInputDataType& rInput)
{
    noalias(rOutput) += rInput;
}

// template instantiations
template void FluidCalculationUtilities::AssignValue<double>(double&, const double&);
template void FluidCalculationUtilities::AssignValue<array_1d<double, 3>>(array_1d<double, 3>&, const array_1d<double, 3>&);
template void FluidCalculationUtilities::AssignValue<Vector>(Vector&, const Vector&);
template void FluidCalculationUtilities::AssignValue<Matrix>(Matrix&, const Matrix&);

template void FluidCalculationUtilities::UpdateValue<array_1d<double, 3>>(array_1d<double, 3>&, const array_1d<double, 3>&);
template void FluidCalculationUtilities::UpdateValue<Vector>(Vector&, const Vector&);
template void FluidCalculationUtilities::UpdateValue<Matrix>(Matrix&, const Matrix&);

} // namespace Kratos
