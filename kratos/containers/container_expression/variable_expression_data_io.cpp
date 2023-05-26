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

template<class TDataType, class TExpressionType, std::size_t... TIndex>
void Assign(
    const TExpressionType& rExpression,
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    if constexpr(
        std::is_same_v<TExpressionType, LiteralFlatExpression<char>> ||
        std::is_same_v<TExpressionType, LiteralFlatExpression<int>>) {
        const auto p_itr = rExpression.cbegin() + EntityDataBeginIndex;
        (rValue.insert_element(TIndex, static_cast<double>(*(p_itr + TIndex))), ...);
    } else if constexpr(std::is_same_v<TExpressionType, LiteralFlatExpression<double>>) {
        const auto p_itr = rExpression.cbegin() + EntityDataBeginIndex;
        (rValue.insert_element(TIndex, *(p_itr + TIndex)), ...);
    } else {
        (rValue.insert_element(TIndex, rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, TIndex)), ...);
    }

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
            << rShape << " ].\n";
    } else if constexpr(std::is_same_v<TDataType, double>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 0)
            << "Constructing a double VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].\n";
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1)
            << "Constructing a vector VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].\n";
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 2)
            << "Constructing a matrix VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].\n";
    } else { // this is for all array_1d types
        KRATOS_ERROR_IF_NOT(rShape == std::vector<IndexType>{std::tuple_size_v<typename TDataType::array_type>})
            << "Constructing a array VariableExpressionDataIO with invalid "
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
template<class TExpressionType>
void VariableExpressionDataIO<TDataType>::Assign(
    TDataType& rOutput,
    const TExpressionType& rExpression,
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
template<class TExpressionType>
void VariableExpressionDataIO<int>::Assign(
    int& rOutput,
    const TExpressionType& rExpression,
    const IndexType EntityIndex) const
{
    if constexpr(
        std::is_same_v<TExpressionType, LiteralFlatExpression<char>> ||
        std::is_same_v<TExpressionType, LiteralFlatExpression<double>>) {
        rOutput = static_cast<int>(*(rExpression.cbegin() + EntityIndex));
    } else if constexpr(std::is_same_v<TExpressionType, LiteralFlatExpression<int>>) {
        rOutput = *(rExpression.cbegin() + EntityIndex);
    } else {
        rOutput = static_cast<int>(rExpression.Evaluate(EntityIndex, EntityIndex, 0));
    }
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
template<class TExpressionType>
void VariableExpressionDataIO<double>::Assign(
    double& rOutput,
    const TExpressionType& rExpression,
    const IndexType EntityIndex) const
{
    if constexpr(
        std::is_same_v<TExpressionType, LiteralFlatExpression<char>> ||
        std::is_same_v<TExpressionType, LiteralFlatExpression<int>>) {
        rOutput = static_cast<double>(*(rExpression.cbegin() + EntityIndex));
    } else if constexpr(std::is_same_v<TExpressionType, LiteralFlatExpression<double>>) {
        rOutput = *(rExpression.cbegin() + EntityIndex);
    } else {
        rOutput = rExpression.Evaluate(EntityIndex, EntityIndex, 0);
    }
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
template<class TExpressionType>
void VariableExpressionDataIO<Vector>::Assign(
    Vector& rOutput,
    const TExpressionType& rExpression,
    const IndexType EntityIndex) const
{
    const IndexType flattened_size = rExpression.GetItemComponentCount();

    if (rOutput.size() != flattened_size) {
        rOutput.resize(flattened_size, false);
    }

    if constexpr(
        std::is_same_v<TExpressionType, LiteralFlatExpression<char>> ||
        std::is_same_v<TExpressionType, LiteralFlatExpression<int>>) {
        const auto p_itr = rExpression.cbegin() + EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput[i] = static_cast<double>(*(p_itr + i));
        }
    } else if constexpr(std::is_same_v<TExpressionType, LiteralFlatExpression<double>>) {
        const auto p_itr = rExpression.cbegin() + EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput[i] = *(p_itr + i);
        }
    } else {
        const IndexType entity_data_begin_index = EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput[i] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
        }
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
template<class TExpressionType>
void VariableExpressionDataIO<Matrix>::Assign(
    Matrix& rOutput,
    const TExpressionType& rExpression,
    const IndexType EntityIndex) const
{
    const auto& r_shape = rExpression.GetItemShape();

    if (rOutput.size1() != r_shape[0] || rOutput.size2() != r_shape[1]) {
        rOutput.resize(r_shape[0], r_shape[1], false);
    }

    const IndexType flattened_size = rExpression.GetItemComponentCount();

    if constexpr(
        std::is_same_v<TExpressionType, LiteralFlatExpression<char>> ||
        std::is_same_v<TExpressionType, LiteralFlatExpression<int>>) {
        const auto p_itr = rExpression.cbegin() + EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput.data()[i] = static_cast<double>(*(p_itr + i));
        }
    } else if constexpr(std::is_same_v<TExpressionType, LiteralFlatExpression<double>>) {
        const auto p_itr = rExpression.cbegin() + EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput.data()[i] = *(p_itr + i);
        }
    } else {
        const IndexType entity_data_begin_index = EntityIndex * flattened_size;
        for (IndexType i = 0; i < flattened_size; ++i) {
            rOutput.data()[i] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
        }
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
#define KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(...) \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionDataIO<__VA_ARGS__>::Assign<LiteralFlatExpression<char>>(__VA_ARGS__&, const LiteralFlatExpression<char>&, const IndexType) const;     \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionDataIO<__VA_ARGS__>::Assign<LiteralFlatExpression<int>>(__VA_ARGS__&, const LiteralFlatExpression<int>&, const IndexType) const;       \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionDataIO<__VA_ARGS__>::Assign<LiteralFlatExpression<double>>(__VA_ARGS__&, const LiteralFlatExpression<double>&, const IndexType) const; \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionDataIO<__VA_ARGS__>::Assign<Expression>(__VA_ARGS__&, const Expression&, const IndexType) const;                                       \
    template class VariableExpressionDataIO<__VA_ARGS__>;

KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(int)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(double)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(array_1d<double, 3>)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(array_1d<double, 4>)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(array_1d<double, 6>)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(array_1d<double, 9>)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(Vector)
KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION(Matrix)

#undef KRATOS_VARIABLE_EXPRESSION_DATA_IO_INSTANTIATION

} // namespace Kratos