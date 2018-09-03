//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_AMATRIX_INTERFACE_H_INCLUDED )
#define  KRATOS_AMATRIX_INTERFACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include "amatrix.h"

// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Internals {

// This class is a modified copy of the AMatrix::Matrix to cover the backward compatibility with ublas in current Kratos implementation.
// The idea is to go step by step toward the original Matrix class and eliminate the differences. So please do not add new API to it.
template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
class Matrix : public AMatrix::MatrixExpression<Matrix<TDataType, TSize1, TSize2>,
                   AMatrix::row_major_access>,
               public AMatrix::MatrixStorage<TDataType, TSize1, TSize2> {
   public:
    using data_type = TDataType;
    using base_type = AMatrix::MatrixStorage<TDataType, TSize1, TSize2>;
    using base_type::at;
    using base_type::data;
    using base_type::size;
    using base_type::size1;
    using base_type::size2;
    using base_type::operator();
    
    using iterator = AMatrix::RandomAccessIterator<TDataType>;
    using const_iterator = AMatrix::RandomAccessIterator<const TDataType>;

    //ublas compatibility definitions
    using value_type = TDataType;
    using size_type = std::size_t;
    using difference_type = std::size_t;
    using const_reference = const TDataType&;
    using reference = TDataType&;
    using const_pointer = TDataType*;
    using pointer = TDataType*;


    Matrix() {}

    explicit Matrix(std::size_t TheSize1, std::size_t TheSize2)
        : base_type(TheSize1, TheSize2) {}

    Matrix(Matrix const& Other) : base_type(Other) {}

    Matrix(Matrix&& Other) : base_type(Other) {}

    template <typename TExpressionType, std::size_t TCategory>
    explicit Matrix(AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other)
        : base_type(Other) {}

    template <typename TOtherMatrixType>
    explicit Matrix(TOtherMatrixType const& Other) : base_type(Other) {}

    explicit Matrix(std::initializer_list<TDataType> InitialValues)
        : base_type(InitialValues) {}

    template <typename TOtherMatrixType>
    Matrix& operator=(TOtherMatrixType const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    Matrix& operator=(Matrix const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    Matrix& operator=(Matrix&& Other) {
        base_type::operator=(Other);
        return *this;
    }

    TDataType& operator()(std::size_t i) { return at(i); }

    TDataType const& operator()(std::size_t i) const { return at(i); }

    friend bool operator==(Matrix const& First, Matrix const& Second) {
        for (std::size_t i = 0; i < First.size(); i++)
            if (First._data[i] != Second._data[i])
                return false;
        return true;
    }

    template <typename TExpressionType, std::size_t TCategory>
    Matrix& operator+=(
        AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other) {
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                at(i, j) += Other(i, j);

        return *this;
    }

    template <typename TExpressionType>
    Matrix& operator+=(
        AMatrix::MatrixExpression<TExpressionType, AMatrix::row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) += Other.expression()[i];

        return *this;
    }

    template <typename TExpressionType, std::size_t TCategory>
    Matrix& operator-=(
        AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other) {
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                at(i, j) -= Other(i, j);

        return *this;
    }

    template <typename TExpressionType>
    Matrix& operator-=(
        AMatrix::MatrixExpression<TExpressionType, AMatrix::row_major_access> const& Other) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) -= Other.expression()[i];

        return *this;
    }

    Matrix& operator*=(data_type TheValue) {
        for (std::size_t i = 0; i < size(); i++)
            at(i) *= TheValue;

        return *this;
    }

    Matrix& operator/=(data_type TheValue) {
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

    iterator begin() { return iterator(data()); }

    iterator end() { return iterator(data() + size()); }

    const_iterator begin() const { return const_iterator(data()); }

    const_iterator end() const { return const_iterator(data() + size()); }

    template <typename TExpressionType>
    data_type dot(
        AMatrix::MatrixExpression<TExpressionType, AMatrix::row_major_access> const& Other)
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

    Matrix& noalias() { return *this; }

    AMatrix::TransposeMatrix<Matrix<TDataType, TSize1, TSize2>> transpose() {
        return AMatrix::TransposeMatrix<Matrix<TDataType, TSize1, TSize2>>(*this);
    }
};

template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
bool operator!=(Matrix<TDataType, TSize1, TSize2> const& First,
    Matrix<TDataType, TSize1, TSize2> const& Second) {
    return !(First == Second);
}

/// output stream function
template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
inline std::ostream& operator<<(std::ostream& rOStream,
    Matrix<TDataType, TSize1, TSize2> const& TheMatrix) {
    rOStream << "{";
    for (std::size_t i = 0; i < TheMatrix.size1(); i++) {
        for (std::size_t j = 0; j < TheMatrix.size2(); j++)
            rOStream << TheMatrix(i, j) << ",";
        rOStream << std::endl;
    }
    rOStream << "}";

    return rOStream;
}

} // namespace Internals

///@name Type Definitions
///@{

using Matrix = Internals::Matrix<double,AMatrix::dynamic, AMatrix::dynamic>;

using Vector = Internals::Matrix<double,AMatrix::dynamic, 1>;

template <typename TDataType> using DenseMatrix=Internals::Matrix<TDataType,AMatrix::dynamic, AMatrix::dynamic>;

template <typename TDataType> using DenseVector=Internals::Matrix<TDataType,AMatrix::dynamic, 1>;

template <typename TDataType, std::size_t TSize> using BoundedVector=Internals::Matrix<TDataType,TSize, 1>;

template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix=Internals::Matrix<TDataType,TSize1, TSize2>;

template <typename TDataType, std::size_t TSize> using BoundedVector=Internals::Matrix<TDataType,TSize, 1>;

template <typename T> T& noalias(T& TheMatrix){return TheMatrix.noalias();}

template <typename T> AMatrix::TransposeMatrix<T> trans(T& TheMatrix){return TheMatrix.transpose();}

using ZeroMatrix = AMatrix::ZeroMatrix<double>;

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
AMatrix::MatrixProductExpression<TExpression1Type, TExpression2Type> prod(
    AMatrix::MatrixExpression<TExpression1Type, TCategory1> const& First,
    AMatrix::MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return AMatrix::MatrixProductExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_AMATRIX_INTERFACE_H_INCLUDED  defined 



