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

template<class TDataType, std::size_t... TIndex>
void ReadVector(
    LiteralFlatExpression<std::conditional_t<std::is_same_v<typename TDataType::value_type, int>, int, double>>& rExpression,
    const IndexType EntityDataBeginIndex,
    const TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    IndexType local_index = 0;
    for (IndexType i = 0; i < rValue.size(); ++i) {
        (rExpression.SetData(EntityDataBeginIndex, local_index++, rValue[i][TIndex]), ...);
    }
}

template<class TDataType, std::size_t... TIndex>
void AssignVector(
    const Expression& rExpression,
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    TDataType& rValue,
    std::index_sequence<TIndex...>)
{
    const auto& r_shape = rExpression.GetItemShape();

    if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                 std::is_same_v<TDataType, array_1d<double, 4>> ||
                 std::is_same_v<TDataType, array_1d<double, 6>> ||
                 std::is_same_v<TDataType, array_1d<double, 9>>) {
        KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 1)
            << "Invalid shape provided. Required shape should have only one dimension. "
            << "Expression dimensions = " << r_shape.size() << ".";
    } else if constexpr(std::is_same_v<TDataType, std::vector<array_1d<double, 3>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 4>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 6>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 9>>>) {
        KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 2)
            << "Invalid shape provided. Required shape should have only 2 dimensions. "
            << "Expression dimensions = " << r_shape.size() << ".";
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }

    if (rValue.size() != r_shape[0]) {
        rValue.resize(r_shape[0]);
    }

    IndexType local_index = 0;
    for (auto& r_value : rValue) {
        (r_value.insert_element(TIndex, rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, local_index++)), ...);
    }
}

} // namespace VariableExpressionDataIOHelperUtilities

template<class TDataType>
VariableExpressionDataIO<TDataType>::VariableExpressionDataIO(const TDataType& SampleValue)
{
    if constexpr(std::is_same_v<TDataType, int> || std::is_same_v<TDataType, double>) {
        mShape.clear();
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        mShape.push_back(SampleValue.size());
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        mShape.push_back(SampleValue.size1());
        mShape.push_back(SampleValue.size2());
    } else if constexpr(
                std::is_same_v<TDataType, array_1d<double, 3>> ||
                std::is_same_v<TDataType, array_1d<double, 4>> ||
                std::is_same_v<TDataType, array_1d<double, 6>> ||
                std::is_same_v<TDataType, array_1d<double, 9>>) {
        mShape.push_back(std::tuple_size_v<typename TDataType::array_type>);
    } else if constexpr(std::is_same_v<TDataType, std::vector<int>>) {
        mShape.push_back(SampleValue.size());
    } else if constexpr(std::is_same_v<TDataType, std::vector<double>>) {
        mShape.push_back(SampleValue.size());
    } else if constexpr(std::is_same_v<TDataType, std::vector<Vector>>) {
        if (SampleValue.size() > 0) {
            mShape.push_back(SampleValue.size());
            mShape.push_back(SampleValue.front().size());
        } else {
            mShape.clear();
        }
    } else if constexpr(std::is_same_v<TDataType, std::vector<Matrix>>) {
        if (SampleValue.size() > 0) {
            mShape.push_back(SampleValue.size());
            mShape.push_back(SampleValue.front().size1());
            mShape.push_back(SampleValue.front().size2());
        } else {
            mShape.clear();
        }
    } else if constexpr(
                std::is_same_v<TDataType, std::vector<array_1d<double, 3>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 4>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 6>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 9>>>) {
        mShape.push_back(SampleValue.size());
        mShape.push_back(std::tuple_size_v<typename TDataType::value_type::array_type>);
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
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
    } else if constexpr(
                std::is_same_v<TDataType, array_1d<double, 3>> ||
                std::is_same_v<TDataType, array_1d<double, 4>> ||
                std::is_same_v<TDataType, array_1d<double, 6>> ||
                std::is_same_v<TDataType, array_1d<double, 9>>) {
        KRATOS_ERROR_IF_NOT(rShape == std::vector<IndexType>{std::tuple_size_v<typename TDataType::array_type>})
            << "Constructing an array VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ] for array of size " << std::tuple_size_v<typename TDataType::array_type> << ".";
    } else if constexpr(std::is_same_v<TDataType, std::vector<int>>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1)
            << "Constructing a vector of ints VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, std::vector<double>>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1)
            << "Constructing a vector of doubles VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, std::vector<Vector>>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 2)
            << "Constructing a vector of vectors VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(std::is_same_v<TDataType, std::vector<Matrix>>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 3)
            << "Constructing a vector of matrices VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ].";
    } else if constexpr(
                std::is_same_v<TDataType, std::vector<array_1d<double, 3>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 4>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 6>>> ||
                std::is_same_v<TDataType, std::vector<array_1d<double, 9>>>) {
        KRATOS_ERROR_IF_NOT(rShape.size() == 2)
            << "Constructing a vector of arrays VariableExpressionDataIO with invalud "
            << "shape. [ shape = "
            << rShape << " ] for array of size " << std::tuple_size_v<typename TDataType::value_type::array_type> << ".";
        KRATOS_ERROR_IF_NOT(rShape[1] == std::tuple_size_v<typename TDataType::value_type::array_type>)
            << "Constructing a vector of arrays VariableExpressionDataIO with invalid "
               "shape. [ shape = "
            << rShape << " ] for array of size " << std::tuple_size_v<typename TDataType::value_type::array_type> << ".";
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
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
    if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                 std::is_same_v<TDataType, array_1d<double, 4>> ||
                 std::is_same_v<TDataType, array_1d<double, 6>> ||
                 std::is_same_v<TDataType, array_1d<double, 9>>) {
        constexpr IndexType N = std::tuple_size_v<typename TDataType::array_type>;
        VariableExpressionDataIOHelperUtilities::Assign(rExpression, EntityIndex, EntityIndex * N, rOutput, std::make_index_sequence<N>{});
    } else if constexpr(std::is_same_v<TDataType, std::vector<array_1d<double, 3>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 4>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 6>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 9>>>) {
        constexpr IndexType N = std::tuple_size_v<typename TDataType::value_type::array_type>;
        VariableExpressionDataIOHelperUtilities::AssignVector(rExpression, EntityIndex, EntityIndex * N * rExpression.GetItemShape()[0], rOutput, std::make_index_sequence<N>{});
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported type.");
    }
}

template<class TDataType>
void VariableExpressionDataIO<TDataType>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const TDataType& Value) const
{
    if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                 std::is_same_v<TDataType, array_1d<double, 4>> ||
                 std::is_same_v<TDataType, array_1d<double, 6>> ||
                 std::is_same_v<TDataType, array_1d<double, 9>>) {
        constexpr IndexType N = std::tuple_size_v<typename TDataType::array_type>;
        VariableExpressionDataIOHelperUtilities::Read(rExpression, EntityIndex * N, Value, std::make_index_sequence<N>{});
    } else if constexpr(std::is_same_v<TDataType, std::vector<array_1d<double, 3>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 4>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 6>>> ||
                        std::is_same_v<TDataType, std::vector<array_1d<double, 9>>>) {
        constexpr IndexType N = std::tuple_size_v<typename TDataType::value_type::array_type>;
        VariableExpressionDataIOHelperUtilities::ReadVector(rExpression, EntityIndex * N * Value.size(), Value, std::make_index_sequence<N>{});
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported type.");
    }
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
void VariableExpressionDataIO<std::vector<int>>::Assign(
    std::vector<int>& rOutput,
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
void VariableExpressionDataIO<int>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const int& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value);
}

template<>
void VariableExpressionDataIO<std::vector<int>>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const std::vector<int>& Value) const
{
    const IndexType entity_data_begin_index = EntityIndex * Value.size();
    for (IndexType i = 0; i < Value.size(); ++i) {
        rExpression.SetData(entity_data_begin_index, i, Value[i]);
    }
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
void VariableExpressionDataIO<std::vector<double>>::Assign(
    std::vector<double>& rOutput,
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
void VariableExpressionDataIO<double>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const double& Value) const
{
    rExpression.SetData(EntityIndex, 0, Value);
}

template<>
void VariableExpressionDataIO<std::vector<double>>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const std::vector<double>& Value) const
{
    const IndexType entity_data_begin_index = EntityIndex * Value.size();
    for (IndexType i = 0; i < Value.size(); ++i) {
        rExpression.SetData(entity_data_begin_index, i, Value[i]);
    }
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
void VariableExpressionDataIO<std::vector<Vector>>::Assign(
    std::vector<Vector>& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    const auto& r_shape = rExpression.GetItemShape();

    KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 2)
        << "Invalid expression size. Requires shape with 2 dimensions. "
        << "Expression dimensions = " << r_shape.size() << ".";

    if (rOutput.size() != r_shape[0]) {
        rOutput.resize(r_shape[0], Vector(r_shape[1]));
    }

    const IndexType flattened_size = r_shape[0] * r_shape[1];
    const IndexType entity_data_begin_index = EntityIndex * flattened_size;

    IndexType local_index = 0;
    for (IndexType i = 0; i < r_shape[0]; ++i) {
        for (IndexType j = 0; j < r_shape[1]; ++j) {
            rOutput[i][j] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, local_index++);
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
void VariableExpressionDataIO<std::vector<Vector>>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const std::vector<Vector>& Value) const
{
    const auto& r_shape = rExpression.GetItemShape();

    KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 2)
        << "Invalid expression size. Requires shape with 2 dimensions. "
        << "Expression dimensions = " << r_shape.size() << ".";

    const IndexType entity_data_begin_index = EntityIndex * r_shape[0] * r_shape[1];
    IndexType local_index = 0;
    for (IndexType i = 0; i < r_shape[0]; ++i) {
        for (IndexType j = 0; j < r_shape[1]; ++j) {
            rExpression.SetData(entity_data_begin_index, local_index++, Value[i][j]);
        }
    }
}

template<>
void VariableExpressionDataIO<Matrix>::Assign(
    Matrix& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    const auto& r_shape = rExpression.GetItemShape();

    KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 2)
        << "Invalid expression size. Requires shape with 2 dimensions. "
        << "Expression dimensions = " << r_shape.size() << ".";

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
void VariableExpressionDataIO<std::vector<Matrix>>::Assign(
    std::vector<Matrix>& rOutput,
    const Expression& rExpression,
    const IndexType EntityIndex) const
{
    const auto& r_shape = rExpression.GetItemShape();

    KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 3)
        << "Invalid expression size. Requires shape with 3 dimensions. "
        << "Expression dimensions = " << r_shape.size() << ".";

    if (rOutput.size() != r_shape[0]) {
        rOutput.resize(r_shape[0], Matrix(r_shape[1], r_shape[2]));
    }

    const IndexType matrix_size = r_shape[1] * r_shape[2];
    const IndexType entity_data_begin_index = EntityIndex * r_shape[0] * matrix_size;

    IndexType local_index = 0;
    for (IndexType i_matrix = 0; i_matrix < r_shape[0]; ++i_matrix) {
        for (IndexType i_component = 0; i_component < matrix_size; ++i_component) {
            rOutput[i_matrix].data()[i_component] = rExpression.Evaluate(EntityIndex, entity_data_begin_index, local_index++);
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

template<>
void VariableExpressionDataIO<std::vector<Matrix>>::Read(
    RawLiteralFlatExpression& rExpression,
    const IndexType EntityIndex,
    const std::vector<Matrix>& Value) const
{
    const auto& r_shape = rExpression.GetItemShape();

    KRATOS_DEBUG_ERROR_IF_NOT(r_shape.size() == 3)
        << "Invalid expression size. Requires shape with 3 dimensions. "
        << "Expression dimensions = " << r_shape.size() << ".";

    const IndexType matrix_size = r_shape[1] * r_shape[2];
    const IndexType entity_data_begin_index = EntityIndex * r_shape[0] * matrix_size;

    IndexType local_index = 0;
    for (IndexType i_matrix = 0; i_matrix < r_shape[0]; ++i_matrix) {
        for (IndexType i_component = 0; i_component < matrix_size; ++i_component) {
            rExpression.SetData(entity_data_begin_index, local_index++, Value[i_matrix].data()[i_component]);
        }
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

template class VariableExpressionDataIO<std::vector<int>>;
template class VariableExpressionDataIO<std::vector<double>>;
template class VariableExpressionDataIO<std::vector<array_1d<double, 3>>>;
template class VariableExpressionDataIO<std::vector<array_1d<double, 4>>>;
template class VariableExpressionDataIO<std::vector<array_1d<double, 6>>>;
template class VariableExpressionDataIO<std::vector<array_1d<double, 9>>>;
template class VariableExpressionDataIO<std::vector<Vector>>;
template class VariableExpressionDataIO<std::vector<Matrix>>;

} // namespace Kratos