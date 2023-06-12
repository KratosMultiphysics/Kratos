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
#include <vector>
#include <utility>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

// Include base h
#include "variable_expression_data_io.h"

namespace Kratos {

namespace VariableExpressionDataIOHelperUtilities {
using IndexType = std::size_t;

template<class TDataType, std::size_t... TIndex>
void Read(
    LiteralFlatExpression<std::conditional_t<std::is_same_v<TDataType, int>, int, double>>& rExpression,
    const IndexType EntityDataBeginIndex,
    const TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    (rExpression.SetData(EntityDataBeginIndex, TIndex, rValue[TIndex]), ...);
}

template<class TDataType, std::size_t... TIndex>
void Assign(
    const Expression& rExpression,
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    (rValue.insert_element(TIndex, rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, TIndex)), ...);
}

} // namespace VariableExpressionDataIOHelperUtilities

template<class TDataType>
VariableExpressionDataIO<TDataType>::VariableExpressionDataIO(const TDataType& SampleValue)
{
    if constexpr(std::is_same_v<TDataType, int>) {
        mShape.clear();
    } else if constexpr(std::is_same_v<TDataType, double>) {
        mShape.clear();
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        mShape.push_back(SampleValue.size());
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        mShape.push_back(SampleValue.size1());
        mShape.push_back(SampleValue.size2());
    } else { // this is for all array_1d types
        mShape.push_back(std::tuple_size_v<typename TDataType::array_type>);
    }
}

template<class TDataType>
VariableExpressionDataIO<TDataType>::VariableExpressionDataIO(const std::vector<IndexType>& rShape)
{
    if constexpr(std::is_same_v<TDataType, int>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 0)
            << "Constructing a int VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, double>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 0)
            << "Constructing a double VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1)
            << "Constructing a vector VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 2)
            << "Constructing a matrix VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else { // this is for all array_1d types
        KRATOS_ERROR_IF_NOT(rShape == std::vector<IndexType>{std::tuple_size_v<typename TDataType::array_type>})
            << "Constructing a array VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ] for array of size " << std::tuple_size_v<typename TDataType::array_type> << ".";
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
    VariableExpressionDataIOHelperUtilities::Assign(rExpression, EntityIndex, EntityIndex * N, rOutput, std::make_index_sequence<N>{});
}

template<class TDataType>
void VariableExpressionDataIO<TDataType>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const TDataType& Value) const
{
    constexpr IndexType N = std::tuple_size_v<typename TDataType::array_type>;
    VariableExpressionDataIOHelperUtilities::Read(rExpression, EntityIndex * N, Value, std::make_index_sequence<N>{});
}

template<>
void VariableExpressionDataIO<int>::Assign(
    int& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    rOutput = static_cast<int>(rExpression.Evaluate(EntityIndex, EntityIndex, 0));
}

template<>
void VariableExpressionDataIO<int>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const int& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value);
}

template<>
void VariableExpressionDataIO<double>::Assign(
    double& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    rOutput = rExpression.Evaluate(EntityIndex, EntityIndex, 0);
}

template<>
void VariableExpressionDataIO<double>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const double& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value);
}

template<>
void VariableExpressionDataIO<Vector>::Assign(
    Vector& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    const IndexType flattened_size = rExpression.GetItemComponentCount();

    if (rOutput.size() != flattened_size) {
        rOutput.resize(flattened_size, false);
    }

    const IndexType entity_data_begin_index = EntityIndex * flattened_size;
    for (IndexType i = 0; i < flattened_size; ++i) {
        rOutput[i] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
    }
}

template<>
void VariableExpressionDataIO<Vector>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const Vector& Value) const
{
    const IndexType entity_data_begin_index = EntityIndex * Value.size();
    for (IndexType i = 0; i < Value.size(); ++i) {
        rExpression.SetData(entity_data_begin_index, i, Value[i]);
    }
}

template<>
void VariableExpressionDataIO<Matrix>::Assign(
    Matrix& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    const auto& r_shape = rExpression.GetItemShape();

    if (rOutput.size1() != r_shape[0] || rOutput.size2() != r_shape[1]) {
        rOutput.resize(r_shape[0], r_shape[1], false);
    }

    const IndexType flattened_size = rExpression.GetItemComponentCount();
    const IndexType entity_data_begin_index = EntityIndex * flattened_size;

    for (IndexType i = 0; i < flattened_size; ++i) {
        rOutput.data()[i] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
    }
}

template<>
void VariableExpressionDataIO<Matrix>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const Matrix& Value) const
{
    const IndexType flattened_size = rExpression.GetItemComponentCount();
    const IndexType entity_data_begin_index = EntityIndex * flattened_size;

    for (IndexType i = 0; i < flattened_size; ++i) {
        rExpression.SetData(entity_data_begin_index, i, Value.data()[i]);
    }
}

// template instantiations
template class VariableExpressionDataIO<int>;
template class VariableExpressionDataIO<double>;
template class VariableExpressionDataIO<array_1d<double, 3>>;
template class VariableExpressionDataIO<array_1d<double, 4>>;
template class VariableExpressionDataIO<array_1d<double, 6>>;
template class VariableExpressionDataIO<array_1d<double, 9>>;
template class VariableExpressionDataIO<Vector>;
template class VariableExpressionDataIO<Matrix>;

} // namespace Kratos