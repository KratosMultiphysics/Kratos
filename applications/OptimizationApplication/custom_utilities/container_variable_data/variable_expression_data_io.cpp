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
#include <utility>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes

// Include base h
#include "variable_expression_data_io.h"

namespace Kratos {

namespace VariableExpressionDataIOHelperUtilities {
using IndexType = std::size_t;

template<class TDataType, std::size_t... TIndex>
void Read(
    LiteralFlatExpression& rExpression,
    const IndexType EntityDataBeginIndex,
    const TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    (rExpression.SetData(EntityDataBeginIndex, TIndex, rValue[TIndex]), ...);
}

template<class TDataType, std::size_t... TIndex>
void Assign(
    const Expression& rExpression,
    const IndexType EntityDataBeginIndex,
    TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    (rValue.insert_element(TIndex, rExpression.Evaluate(EntityDataBeginIndex, TIndex)), ...);
}

} // namespace VariableExpressionDataIOHelperUtilities

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

template<class TDataType>
std::shared_ptr<VariableExpressionDataIO<TDataType>> VariableExpressionDataIO<TDataType>::Create(const TDataType& SampleValue)
{
    return std::make_shared<VariableExpressionDataIO<TDataType>>(SampleValue);
}

template<class TDataType>
std::shared_ptr<VariableExpressionDataIO<TDataType>> VariableExpressionDataIO<TDataType>::Create(const std::vector<IndexType>& rShape)
{
    return std::make_shared<VariableExpressionDataIO<TDataType>>(rShape);
}

template<class TDataType>
void VariableExpressionDataIO<TDataType>::Assign(
    TDataType& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    constexpr IndexType N = std::tuple_size_v<typename TDataType::array_type>;
    VariableExpressionDataIOHelperUtilities::Assign(rExpression, EntityIndex * N, rOutput, std::make_index_sequence<N>{});
}

template<class TDataType>
void VariableExpressionDataIO<TDataType>::Read(
    LiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const TDataType& Value) const
{
    constexpr IndexType N = std::tuple_size_v<typename TDataType::array_type>;
    VariableExpressionDataIOHelperUtilities::Read(rExpression, EntityIndex * N, Value, std::make_index_sequence<N>{});
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

// template instantiations
template class VariableExpressionDataIO<double>;
template class VariableExpressionDataIO<array_1d<double, 3>>;

} // namespace Kratos