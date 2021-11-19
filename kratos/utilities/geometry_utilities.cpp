//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes

// Include base h
#include "geometry_utilities.h"

namespace Kratos
{
template <class TDataType>
void GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint(
    TDataType& rOutput,
    const GeometryType& rGeometry,
    const Variable<TDataType>& rVariable,
    const Vector& rGaussPointShapeFunctionValues,
    const int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = rGeometry.PointsNumber();

    noalias(rOutput) = rGeometry[0].FastGetSolutionStepValue(rVariable, Step) *
                       rGaussPointShapeFunctionValues[0];

    for (SizeType i_node = 1; i_node < number_of_nodes; ++i_node)
    {
        noalias(rOutput) += rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step) *
                            rGaussPointShapeFunctionValues[i_node];
    }

    KRATOS_CATCH("");
}

template <>
void KRATOS_API(KRATOS_CORE) GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint<double>(
    double& rOutput,
    const GeometryType& rGeometry,
    const Variable<double>& rVariable,
    const Vector& rGaussPointShapeFunctionValues,
    const int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = rGeometry.PointsNumber();

    rOutput = rGeometry[0].FastGetSolutionStepValue(rVariable, Step) *
              rGaussPointShapeFunctionValues[0];
    for (SizeType i_node = 1; i_node < number_of_nodes; ++i_node)
    {
        rOutput += rGeometry[i_node].FastGetSolutionStepValue(rVariable, Step) *
                   rGaussPointShapeFunctionValues[i_node];
    }

    KRATOS_CATCH("");
}

void GeometryUtils::EvaluateHistoricalVariableGradientAtGaussPoint(
    array_1d<double, 3>& rOutput,
    const GeometryType& rGeometry,
    const Variable<double>& rVariable,
    const Matrix& rGaussPointShapeFunctionDerivativeValues,
    const int Step)
{
    noalias(rOutput) = ZeroVector(3);
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension = rGaussPointShapeFunctionDerivativeValues.size2();

    for (SizeType a = 0; a < number_of_nodes; ++a)
    {
        const double value = rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (SizeType i = 0; i < dimension; ++i)
            rOutput[i] += rGaussPointShapeFunctionDerivativeValues(a, i) * value;
    }
}

void GeometryUtils::EvaluateHistoricalVariableGradientAtGaussPoint(
    BoundedMatrix<double, 3, 3>& rOutput,
    const GeometryType& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rGaussPointShapeFunctionDerivativeValues,
    const int Step)
{
    noalias(rOutput) = ZeroMatrix(3, 3);
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension = rGaussPointShapeFunctionDerivativeValues.size2();

    for (SizeType a = 0; a < number_of_nodes; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (SizeType i = 0; i < dimension; ++i)
        {
            for (SizeType j = 0; j < dimension; ++j)
            {
                rOutput(i, j) +=
                    rGaussPointShapeFunctionDerivativeValues(a, j) * r_value[i];
            }
        }
    }
}

// template instantiations

template void KRATOS_API(KRATOS_CORE) GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint<array_1d<double, 3>>(
    array_1d<double, 3>& rOutput,
    const GeometryType&,
    const Variable<array_1d<double, 3>>&,
    const Vector&,
    const int);

} // namespace Kratos.
