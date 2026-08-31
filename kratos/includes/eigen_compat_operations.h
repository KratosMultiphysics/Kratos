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
#include <boost/numeric/ublas/expression_types.hpp>

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

/// Matrix/matrix and matrix/vector product (lazy Eigen expression). A vector
/// first operand keeps the ublas prod(v, M) semantics (v^T M), returned
/// column-shaped so it assigns to the vector types.
template<class TDerived1, class TDerived2>
inline auto prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    if constexpr (TDerived1::ColsAtCompileTime == 1 &&
                  TDerived2::ColsAtCompileTime != 1 &&
                  TDerived2::RowsAtCompileTime != 1) { // a row-shaped rB is a plain (outer) product
        return (rA.transpose() * rB).transpose();
    } else {
        return rA * rB;
    }
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

/// uBLAS defines trans() on a *vector* as the identity ((trans v)[i] = v[i],
/// boost vector_expression.hpp), and it discriminates vector from matrix by
/// TYPE, not by shape: an N x 1 matrix transposes, a vector does not. The
/// concrete Kratos vector types keep that semantic here (the bounded vector
/// and array_1d counterparts live in eigen_ublas_compat_operations.h), so
/// idioms like outer_prod(x, trans(y)) and prod(trans(v), M) mean the same
/// under both backends.
template<class TDataType>
inline const EigenVector<TDataType>& trans(const EigenVector<TDataType>& rV)
{
    return rV;
}

/// Transpose of a sparse matrix — e.g. for 0.5 * (K + trans(K)) symmetrization
/// of a system matrix. Unlike the dense overload this MATERIALIZES the result
/// (one O(nnz) copy): Eigen's lazy sparse transpose flips the storage order,
/// and sparse binary operations require both sides to share the same order.
template<class TDerived>
inline typename TDerived::PlainObject trans(const Eigen::SparseMatrixBase<TDerived>& rM)
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

/// Frobenius norm (same as norm_2 for vectors, as in ublas).
template<class TDerived>
inline auto norm_frobenius(const Eigen::MatrixBase<TDerived>& rM)
{
    return rM.norm();
}

/// Element-wise product (lazy Eigen expression).
template<class TDerived1, class TDerived2>
inline auto element_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.cwiseProduct(rY.derived());
}

/// Element-wise division (lazy Eigen expression).
template<class TDerived1, class TDerived2>
inline auto element_div(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.cwiseQuotient(rY.derived());
}

namespace Internals {

/// noalias proxy for the dynamic Eigen-backed types. Eigen RHS go through
/// Eigen's own NoAlias (the exact ublas noalias contract: assignment without
/// a protective temporary); uBLAS expression RHS (ZeroMatrix, IdentityMatrix,
/// prod(...) of uBLAS operands, ...) are evaluated element-wise, resizing the
/// target as plain ublas assignment would.
template<class TTarget>
class EigenDynamicNoAliasProxy
{
public:
    explicit EigenDynamicNoAliasProxy(TTarget& rTarget) : mrTarget(rTarget) {}

    template<class TDerived>
    TTarget& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() = rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TTarget& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() += rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TTarget& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() -= rExpression.derived();
        return mrTarget;
    }

    template<class TExpression>
    TTarget& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        mrTarget.resize(r_e.size1(), r_e.size2(), false);
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                mrTarget(i, j) = r_e(i, j);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                mrTarget(i, j) += r_e(i, j);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                mrTarget(i, j) -= r_e(i, j);
        return mrTarget;
    }

    template<class TExpression>
    TTarget& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        mrTarget.resize(r_e.size(), false);
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] = r_e(i);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] += r_e(i);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] -= r_e(i);
        return mrTarget;
    }

private:
    TTarget& mrTarget;
};

} // namespace Internals

/// noalias proxies accepting both Eigen and uBLAS right-hand sides.
template<class TDataType>
inline auto noalias(EigenMatrix<TDataType>& rM)
{
    return Internals::EigenDynamicNoAliasProxy<EigenMatrix<TDataType>>(rM);
}

template<class TDataType>
inline auto noalias(EigenVector<TDataType>& rV)
{
    return Internals::EigenDynamicNoAliasProxy<EigenVector<TDataType>>(rV);
}

/// Row proxy (readable and writable, as ublas row()).
/// The uBLAS matrix_row models a *vector*; Eigen's vector convention is a
/// column, so the row block is transposed to keep vector semantics
/// (assignment to array_1d/EigenVector, inner_prod, v -= row(M,i), ...).
template<class TDataType>
inline auto row(EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
}

template<class TDataType>
inline auto row(const EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
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

/// Matrix subrange [Row1, Row2) x [Col1, Col2) block (readable and writable, as ublas subrange()).
template<class TDataType>
inline auto subrange(EigenMatrix<TDataType>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

template<class TDataType>
inline auto subrange(const EigenMatrix<TDataType>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

} // namespace Kratos
