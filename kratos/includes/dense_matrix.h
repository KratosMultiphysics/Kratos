#pragma once

// A matrix Library to be simple and fast
#include <cmath>
#include <limits>
#include "matrix_storage.h"

namespace Kratos {

template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
class DenseMatrix : public AMatrix::MatrixExpression<DenseMatrix<TDataType, TSize1, TSize2>,
                   row_major_access>,
               public MatrixStorage<TDataType, TSize1, TSize2> {
   public:
    using data_type = TDataType;
    using base_type = MatrixStorage<TDataType, TSize1, TSize2>;
    using base_type::at;
    using base_type::size;
    using base_type::size1;
    using base_type::size2;
    using base_type::operator();

    DenseMatrix() {}

    explicit DenseMatrix(std::size_t TheSize1, std::size_t TheSize2)
        : base_type(TheSize1, TheSize2) {}

    DenseMatrix(DenseMatrix const& Other) : base_type(Other) {}

    DenseMatrix(DenseMatrix&& Other) : base_type(Other) {}

    template <typename TExpressionType, std::size_t TCategory>
    explicit DenseMatrix(AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other)
        : base_type(Other) {}

    template <typename TOtherMatrixType>
    explicit DenseMatrix(TOtherMatrixType const& Other) : base_type(Other) {}

    explicit DenseMatrix(std::initializer_list<TDataType> InitialValues)
        : base_type(InitialValues) {}

    template <typename TOtherMatrixType>
    DenseMatrix& operator=(TOtherMatrixType const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    DenseMatrix& operator=(DenseMatrix const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    DenseMatrix& operator=(DenseMatrix&& Other) {
        base_type::operator=(Other);
        return *this;
    }

    TDataType& operator()(std::size_t i) { return at(i); }

    TDataType const& operator()(std::size_t i) const { return at(i); }

    friend bool operator==(DenseMatrix const& First, DenseMatrix const& Second) {
        for (std::size_t i = 0; i < First.size(); i++)
            if (First._data[i] != Second._data[i])
                return false;
        return true;
    }

    template <typename TExpressionType, std::size_t TCategory>
    DenseMatrix& operator+=(
        AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other) {
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                at(i, j) += Other(i, j);

        return *this;
    }

    template <typename TExpressionType>
    DenseMatrix& operator+=(
        AMatrix::MatrixExpression<TExpressionType, row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) += Other.expression()[i];

        return *this;
    }

    template <typename TExpressionType, std::size_t TCategory>
    DenseMatrix& operator-=(
        AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other) {
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                at(i, j) -= Other(i, j);

        return *this;
    }

    template <typename TExpressionType>
    DenseMatrix& operator-=(
        AMatrix::MatrixExpression<TExpressionType, row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) -= Other.expression()[i];

        return *this;
    }

    DenseMatrix& operator*=(data_type TheValue) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) *= TheValue;

        return *this;
    }

    DenseMatrix& operator/=(data_type TheValue) {
        auto inverse_of_value = 1.00 / TheValue;
        for (std::size_t i = 0; i < size(); i++)
            at(i) *= inverse_of_value;

        return *this;
    }


    void resize(std::size_t NewSize1, std::size_t NewSize2, bool preserve = 0){
        KRATOS_DEBUG_ERROR_IF(preserve) << "The preserve is not supported anymore" << std::endl;

        base_type::resize(NewSize1,NewSize2);
    }

    void resize(std::size_t NewSize, bool preserve = 0){
        KRATOS_DEBUG_ERROR_IF(preserve) << "The preserve is not supported anymore" << std::endl;
        base_type::resize(NewSize);
    }

    template <typename TExpressionType>
    data_type dot(
        AMatrix::MatrixExpression<TExpressionType, row_major_access> const& Other)
        const {
        data_type result = data_type();
        for (std::size_t i = 0; i < size(); ++i) {
            result += at(i) * Other.expression()[i];
        }
        return result;
    }

    data_type squared_norm() const { return dot(*this); }

    data_type norm() const { return std::sqrt(dot(*this)); }

    void normalize() {
        auto the_norm = norm();
        if (the_norm > std::numeric_limits<data_type>::epsilon()) {
            const auto norm_inverse = 1.0 / the_norm;
            for (std::size_t i = 0; i < size(); ++i) {
                at(i) *= norm_inverse;
            }
        }
    }

    DenseMatrix& noalias() { return *this; }

    TransposeMatrix<DenseMatrix<TDataType, TSize1, TSize2>> transpose() {
        return TransposeMatrix<DenseMatrix<TDataType, TSize1, TSize2>>(*this);
    }
};

template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
bool operator!=(DenseMatrix<TDataType, TSize1, TSize2> const& First,
    DenseMatrix<TDataType, TSize1, TSize2> const& Second) {
    return !(First == Second);
}

/// output stream function
template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
inline std::ostream& operator<<(std::ostream& rOStream,
    DenseMatrix<TDataType, TSize1, TSize2> const& TheMatrix) {
    rOStream << "{";
    for (std::size_t i = 0; i < TheMatrix.size1(); i++) {
        for (std::size_t j = 0; j < TheMatrix.size2(); j++)
            rOStream << TheMatrix(i, j) << ",";
        rOStream << std::endl;
    }
    rOStream << "}";

    return rOStream;
}

}  // namespace Kratos
