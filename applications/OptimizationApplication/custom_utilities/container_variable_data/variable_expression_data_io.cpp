//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes

// Include base h
#include "variable_expression_data_io.h"

namespace Kratos {

template<class TDataType>
VariableExpressionDataIO<TDataType>::VariableExpressionDataIO(const TDataType& SampleValue)
{
    if constexpr(std::is_same_v<TDataType, double>) {
        mShape.clear();
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        mShape.push_back(3);
    }
}

template<class TDataType>
VariableExpressionDataIO<TDataType>::VariableExpressionDataIO(const std::vector<IndexType>& rShape)
{
    if constexpr(std::is_same_v<TDataType, double>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 0)
            << "Constructing a double VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].\n";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        KRATOS_ERROR_IF_NOT(rShape == std::vector<IndexType>{3})
            << "Constructing a array3 VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].\n";
    }
}

template<>
void VariableExpressionDataIO<double>::Assign(
    double& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    rOutput = rExpression.Evaluate(EntityIndex, 0);
}

template<>
void VariableExpressionDataIO<double>::Read(
    LiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const double& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value);
}

template<>
void VariableExpressionDataIO<array_1d<double, 3>>::Assign(
    array_1d<double, 3>& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    rOutput[0] = rExpression.Evaluate(EntityIndex, 0);
    rOutput[1] = rExpression.Evaluate(EntityIndex, 1);
    rOutput[2] = rExpression.Evaluate(EntityIndex, 2);
}

template<>
void VariableExpressionDataIO<array_1d<double, 3>>::Read(
    LiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const array_1d<double, 3>& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value[0]);
    rExpression.SetData(EntityIndex, 1, Value[1]);
    rExpression.SetData(EntityIndex, 2, Value[2]);
}

// template instantiations
template class VariableExpressionDataIO<double>;
template class VariableExpressionDataIO<array_1d<double, 3>>;

} // namespace Kratos