#pragma once

#include "matrix_expression.h"

namespace AMatrix {

template <typename TDataType, std::size_t TSize>
class DenseStorage {
    TDataType _data[TSize];

   public:
    DenseStorage() {}

    explicit DenseStorage(std::size_t TheSize) {}

    DenseStorage(DenseStorage const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
    }

    DenseStorage(DenseStorage&& Other) = default;

    template <typename TExpressionType, std::size_t TCategory>
    explicit DenseStorage(
        MatrixExpression<TExpressionType, TCategory> const& Other)
        : DenseStorage(Other.expression()) {}

    template <typename TExpressionType>
    explicit DenseStorage(
        MatrixExpression<TExpressionType, row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other.expression()[i];
    }

    template <typename TOtherMatrixType>
    explicit DenseStorage(TOtherMatrixType const& Other) {
        auto i_data = _data;
        for (std::size_t i = 0; i < Other.size1(); i++)
            for (std::size_t j = 0; j < Other.size2(); j++)
                *(i_data++) = Other(i, j);
    }

    explicit DenseStorage(std::initializer_list<TDataType> InitialValues) {
        std::size_t position = 0;
        for (auto& i : InitialValues) {
            _data[position++] = i;
        }
    }

    template <typename TExpressionType, std::size_t TCategory>
    DenseStorage& operator=(
        MatrixExpression<TExpressionType, TCategory> const& Other) {
        auto i_data = _data;
        for (std::size_t i = 0; i <  Other.expression().size1(); i++)
            for (std::size_t j = 0; j <  Other.expression().size2(); j++)
                *(i_data++) = Other.expression()(i, j);
        return *this;
    }

    template <typename TExpressionType>
    DenseStorage& operator=(
        MatrixExpression<TExpressionType, row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other.expression()[i];
        return *this;
    }

    template <typename TOtherMatrixType>
    DenseStorage& operator=(TOtherMatrixType const& Other) {
        auto i_data = _data;
        for (std::size_t i = 0; i < Other.size1(); i++)
            for (std::size_t j = 0; j < Other.size2(); j++)
                *(i_data++) = Other(i, j);
        return *this;
    }

    DenseStorage& operator=(DenseStorage const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
        return *this;
    }

    DenseStorage& operator=(DenseStorage&& Other) = default;

    TDataType& operator[](std::size_t i) { return at(i); }

    TDataType const& operator[](std::size_t i) const { return at(i); }

    TDataType& at(std::size_t i) { return _data[i]; }

    TDataType const& at(std::size_t i) const { return _data[i]; }

    static constexpr std::size_t size() { return TSize; }

    TDataType* data() { return _data; }

    TDataType const* data() const { return _data; }
};

} // namespace AMatrix