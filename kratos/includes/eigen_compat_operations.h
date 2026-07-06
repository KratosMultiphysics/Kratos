//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <cstddef>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes

namespace Kratos {

// Forward declarations of the wrapper types (defined in kratos_eigen_interface.h,
// which includes this header at its end).
template<class TDataType> class EigenMatrix;
template<class TDataType> class EigenVector;

///@name Eigen compatibility operations
///@{
// uBLAS-style free functions for Eigen types, so that generic Kratos code
// written with the ublas idiom (prod, inner_prod, noalias, trans, ...) also
// compiles against the Eigen backend types.
//
// Two overload families are needed to coexist with the boost::numeric::ublas
// templates injected into namespace Kratos by ublas_interface.h:
// - Functions taking *expressions* (prod, inner_prod, trans, norms, ...)
//   deduce on Eigen::MatrixBase / Eigen::SparseMatrixBase; the ublas
//   counterparts deduce on ublas::matrix_expression / vector_expression, so
//   each library's overloads drop out of the candidate set for the other's
//   argument types.
// - Functions taking a *mutable target* (noalias, row, column, subrange)
//   must overload on the concrete Kratos wrapper types: the ublas templates
//   for these are fully generic (e.g. template<class C> noalias(C&)) and
//   would otherwise be an exact match for the Eigen types too. Overloading on
//   EigenMatrix<T>/EigenVector<T> wins by partial ordering.
///@}

/// Matrix/matrix and matrix/vector product (lazy Eigen expression).
template<class TDerived1, class TDerived2>
inline auto prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return rA * rB;
}

/// Sparse-matrix/vector (or sparse/dense) product.
template<class TDerived1, class TDerived2>
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rX)
{
    return rA * rX;
}

/// Sparse-matrix/sparse-matrix product.
template<class TDerived1, class TDerived2>
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::SparseMatrixBase<TDerived2>& rB)
{
    return rA * rB;
}

/// Scalar product of two vectors.
template<class TDerived1, class TDerived2>
inline auto inner_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.dot(rY);
}

/// Outer product of two vectors.
template<class TDerived1, class TDerived2>
inline auto outer_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX * rY.transpose();
}

/// Transpose (lazy Eigen expression).
template<class TDerived>
inline auto trans(const Eigen::MatrixBase<TDerived>& rM)
{
    return rM.transpose();
}

/// 1-norm.
template<class TDerived>
inline auto norm_1(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.template lpNorm<1>();
}

/// Euclidean norm (Frobenius norm for matrices, as in ublas).
template<class TDerived>
inline auto norm_2(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.norm();
}

/// Infinity norm.
template<class TDerived>
inline auto norm_inf(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.template lpNorm<Eigen::Infinity>();
}

/// Sum of all coefficients.
template<class TDerived>
inline auto sum(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.sum();
}

/// noalias proxy: Eigen's own NoAlias supports =, += and -=, which is exactly
/// the ublas noalias contract (assignment without a protective temporary).
template<class TDataType>
inline auto noalias(EigenMatrix<TDataType>& rM)
{
    return rM.noalias();
}

template<class TDataType>
inline auto noalias(EigenVector<TDataType>& rV)
{
    return rV.noalias();
}

/// Row proxy (readable and writable, as ublas row()).
template<class TDataType>
inline auto row(EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I);
}

template<class TDataType>
inline auto row(const EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I);
}

/// Column proxy (readable and writable, as ublas column()).
template<class TDataType>
inline auto column(EigenMatrix<TDataType>& rM, const std::size_t J)
{
    return rM.col(J);
}

template<class TDataType>
inline auto column(const EigenMatrix<TDataType>& rM, const std::size_t J)
{
    return rM.col(J);
}

/// Vector subrange [Low, High) proxy (readable and writable, as ublas subrange()).
template<class TDataType>
inline auto subrange(EigenVector<TDataType>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

template<class TDataType>
inline auto subrange(const EigenVector<TDataType>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

} // namespace Kratos
