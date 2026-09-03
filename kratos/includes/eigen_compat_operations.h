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

/**
 * @brief Eigen compatibility operations 
 * @details uBLAS-style free functions for Eigen types, so that generic Kratos code
 * written with the ublas idiom (prod, inner_prod, noalias, trans, ...) also
 * compiles against the Eigen backend types.
 *
 * Two overload families are needed to coexist with the boost::numeric::ublas
 * templates injected into namespace Kratos by ublas_interface.h:
 * - Functions taking *expressions* (prod, inner_prod, trans, norms, ...)
 *   deduce on Eigen::MatrixBase / Eigen::SparseMatrixBase; the ublas
 *   counterparts deduce on ublas::matrix_expression / vector_expression, so
 *   each library's overloads drop out of the candidate set for the other's
 *   argument types.
 * - Functions taking a *mutable target* (noalias, row, column, subrange)
 *   must overload on the concrete Kratos wrapper types: the ublas templates
 *   for these are fully generic (e.g. template<class C> noalias(C&)) and
 *   would otherwise be an exact match for the Eigen types too. Overloading on
 *   EigenMatrix<T>/EigenVector<T> wins by partial ordering.
 */

namespace Internals {

/**
 * @brief Shared implementation for the dense prod() overloads below.
 * @details Kept out of the "prod" overload set on purpose: a call from
 * inside a prod() overload to unqualified prod(rA, rB) would re-run overload
 * resolution against *every* prod overload (including the sparse ones), and
 * deducing those against a reference-to-MatrixBase argument forces an
 * instantiation of Eigen::SparseMatrixBase<TDerived> that hard-errors (no
 * Eigen::internal::traits specialization exists for the dense Kratos wrapper
 * types).
 * @tparam TDerived1 Concrete type of the left operand.
 * @tparam TDerived2 Concrete type of the right operand.
 * @param rA Left-hand operand.
 * @param rB Right-hand operand.
 * @return The product as a lazy Eigen expression; a row-shaped rB paired
 * with a column-shaped rA is computed as an outer product instead.
 */
template<class TDerived1, class TDerived2>
inline auto dense_prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    if constexpr (TDerived1::ColsAtCompileTime == 1 &&
                  TDerived2::ColsAtCompileTime != 1 &&
                  TDerived2::RowsAtCompileTime != 1) { // a row-shaped rB is a plain (outer) product
        return (rA.transpose() * rB).transpose();
    } else {
        return rA * rB;
    }
}

} // namespace Internals

/**
 * @brief Matrix/matrix and matrix/vector product (lazy Eigen expression).
 * @details A vector first operand keeps the ublas prod(v, M) semantics
 * (v^T M), returned column-shaped so it assigns to the vector types.
 *
 * The requires-clause guards TDerived1/TDerived2 the same way as the sparse
 * overloads below: prod<TResult>(...) shares the "prod" name, so a call
 * like prod<Matrix>(...) explicitly substitutes TDerived1 = Matrix into
 * *every* prod overload, this one included, and naming
 * Eigen::MatrixBase<Matrix> for that substitution forces Eigen to complete
 * it — which hard-errors (Eigen::internal::traits<Matrix> is unspecialized
 * for the dense Kratos wrapper types). A genuinely deduced TDerived1/2 (no
 * explicit argument given) never hits this: deduction against a concrete
 * Kratos::Matrix/Vector/BoundedXxx argument always resolves to its Eigen
 * *base* (Eigen::Matrix<...>, which Eigen itself specializes traits for),
 * never to the Kratos wrapper type itself.
 * @tparam TDerived1 Concrete type of the left operand.
 * @tparam TDerived2 Concrete type of the right operand.
 * @param rA Left-hand operand.
 * @param rB Right-hand operand.
 * @return The product as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
    requires requires {
        typename Eigen::internal::traits<TDerived1>::StorageKind;
        typename Eigen::internal::traits<TDerived2>::StorageKind;
    }
inline auto prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return Internals::dense_prod(rA, rB);
}

/**
 * @brief Matrix/matrix and matrix/vector product with an explicitly pinned result type.
 * @details uBLAS lets prod() pin the result type explicitly (prod<Vector>(A, B),
 * prod<Matrix>(A, B), ...), e.g. to combine a vector-shaped and a
 * matrix-shaped operand in one expression without an intermediate local.
 * TResult only appears in the return type, so it is never deduced and this
 * overload only participates when the caller supplies it explicitly;
 * otherwise the plain two-argument overload above is the sole match.
 * @tparam TResult Result type explicitly requested by the caller.
 * @tparam TDerived1 Concrete type of the left operand.
 * @tparam TDerived2 Concrete type of the right operand.
 * @param rA Left-hand operand.
 * @param rB Right-hand operand.
 * @return The product, evaluated into a TResult.
 */
template<class TResult, class TDerived1, class TDerived2>
inline TResult prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return TResult(Internals::dense_prod(rA, rB));
}

/**
 * @brief Sparse-matrix/vector (or sparse/dense) product.
 * @details The requires-clause guards TDerived1 before it is ever used to
 * name Eigen::SparseMatrixBase<TDerived1>: prod<TResult>(...) (above) shares
 * the "prod" name, so a call like prod<Matrix>(...) explicitly substitutes
 * TDerived1 = Matrix into *every* prod overload, this one included. Naming
 * SparseMatrixBase<Matrix> for that substitution forces Eigen to complete
 * it, which needs Eigen::internal::traits<Matrix> — unspecialized for the
 * dense Kratos wrapper types — and that incompleteness surfaces several
 * instantiations away from this "immediate context", as a hard error rather
 * than a SFINAE-friendly one. Checking the (specialized-or-not) traits
 * directly first, before naming SparseMatrixBase<TDerived1>, keeps the
 * failure a one-step, immediate-context substitution failure instead.
 * @tparam TDerived1 Concrete type of the sparse left operand.
 * @tparam TDerived2 Concrete type of the dense right operand.
 * @param rA Left-hand (sparse) operand.
 * @param rX Right-hand (dense) operand.
 * @return The product as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
    requires requires { typename Eigen::internal::traits<TDerived1>::StorageKind; }
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rX)
{
    return rA * rX;
}

/**
 * @brief Sparse-matrix/sparse-matrix product.
 * @details See the requires-clause guard note on the sparse/dense overload above.
 * @tparam TDerived1 Concrete type of the left operand.
 * @tparam TDerived2 Concrete type of the right operand.
 * @param rA Left-hand operand.
 * @param rB Right-hand operand.
 * @return The product as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
    requires requires {
        typename Eigen::internal::traits<TDerived1>::StorageKind;
        typename Eigen::internal::traits<TDerived2>::StorageKind;
    }
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::SparseMatrixBase<TDerived2>& rB)
{
    return rA * rB;
}

/**
 * @brief uBLAS-style axpy product: writes the product into a target instead of returning an expression.
 * @details Computes rY = A B for Init = true, rY += A B otherwise. A single
 * overload covers the three uBLAS forms, dispatching on the compile-time
 * shape of the first operand exactly as prod() does above: matrix x vector,
 * matrix x matrix, and — for a vector first operand — the vector x matrix
 * form, which keeps the ublas axpy_prod(v, M, y) semantics (y = M^T v).
 *
 * As in uBLAS the target must not alias either operand: the assignment goes
 * through Eigen's noalias() path, so no protective temporary is built. That
 * same path resizes rY when it is a resizable plain object, matching the
 * ublas containers' behaviour for the generic Kratos code that calls
 * axpy_prod (UblasSpace::Mult / ::TransposeMult instantiated on the
 * Eigen-backed Kratos::Matrix/Vector, for instance).
 * @tparam TDerived1 Concrete type of the first operand (the vector, in the vector x matrix form).
 * @tparam TDerived2 Concrete type of the second operand (the matrix, in the vector x matrix form).
 * @tparam TTargetType Concrete type of the target (a resizable Eigen-backed matrix or vector).
 * @param rA First operand.
 * @param rB Second operand.
 * @param rY Target the product is written into.
 * @param Init If true (default) rY is overwritten; if false the product is added to rY.
 * @return Reference to rY.
 */
template<class TDerived1, class TDerived2, class TTargetType>
inline TTargetType& axpy_prod(
    const Eigen::MatrixBase<TDerived1>& rA,
    const Eigen::MatrixBase<TDerived2>& rB,
    TTargetType& rY,
    const bool Init = true)
{
    if constexpr (TDerived1::ColsAtCompileTime == 1) { // vector first operand: axpy_prod(v, M, y) is y = M^T v
        if (Init) {
            rY.noalias() = rB.derived().transpose() * rA.derived();
        } else {
            rY.noalias() += rB.derived().transpose() * rA.derived();
        }
    } else {
        if (Init) {
            rY.noalias() = rA.derived() * rB.derived();
        } else {
            rY.noalias() += rA.derived() * rB.derived();
        }
    }
    return rY;
}

/**
 * @brief Sparse counterpart of the axpy product above (sparse matrix times a dense vector or matrix).
 * @details See the note on the dense overload above regarding aliasing and resizing.
 * @tparam TDerived1 Concrete type of the sparse left operand.
 * @tparam TDerived2 Concrete type of the dense right operand.
 * @tparam TTargetType Concrete type of the target (a resizable Eigen-backed matrix or vector).
 * @param rA Left-hand (sparse) operand.
 * @param rX Right-hand (dense) operand.
 * @param rY Target the product is written into.
 * @param Init If true (default) rY is overwritten; if false the product is added to rY.
 * @return Reference to rY.
 */
template<class TDerived1, class TDerived2, class TTargetType>
inline TTargetType& axpy_prod(
    const Eigen::SparseMatrixBase<TDerived1>& rA,
    const Eigen::MatrixBase<TDerived2>& rX,
    TTargetType& rY,
    const bool Init = true)
{
    if (Init) {
        rY.noalias() = rA.derived() * rX.derived();
    } else {
        rY.noalias() += rA.derived() * rX.derived();
    }
    return rY;
}

/**
 * @brief Scalar (dot) product of two vectors.
 * @tparam TDerived1 Concrete type of the first vector expression.
 * @tparam TDerived2 Concrete type of the second vector expression.
 * @param rX First vector.
 * @param rY Second vector.
 * @return The dot product of rX and rY.
 */
template<class TDerived1, class TDerived2>
inline auto inner_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.dot(rY);
}

/**
 * @brief Outer product of two vectors.
 * @tparam TDerived1 Concrete type of the first (column) vector expression.
 * @tparam TDerived2 Concrete type of the second (column) vector expression.
 * @param rX First vector, forming the column space of the result.
 * @param rY Second vector, transposed to form the row space of the result.
 * @return The outer product rX * rY^T as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
inline auto outer_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX * rY.transpose();
}

/**
 * @brief Transpose of a dense matrix or vector expression.
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rM Expression to transpose.
 * @return The transpose as a lazy Eigen expression.
 */
template<class TDerived>
inline auto trans(const Eigen::MatrixBase<TDerived>& rM)
{
    return rM.transpose();
}

/**
 * @brief Transpose of an EigenVector: the identity, per uBLAS vector semantics.
 * @details uBLAS defines trans() on a *vector* as the identity ((trans v)[i] = v[i],
 * boost vector_expression.hpp), and it discriminates vector from matrix by
 * TYPE, not by shape: an N x 1 matrix transposes, a vector does not. The
 * concrete Kratos vector types keep that semantic here (the bounded vector
 * and array_1d counterparts live in eigen_ublas_compat_operations.h), so
 * idioms like outer_prod(x, trans(y)) and prod(trans(v), M) mean the same
 * under both backends.
 * @tparam TDataType Scalar type of the vector.
 * @param rV Vector to "transpose".
 * @return Reference to rV, unchanged.
 */
template<class TDataType>
inline const EigenVector<TDataType>& trans(const EigenVector<TDataType>& rV)
{
    return rV;
}

/**
 * @brief Transpose of a sparse matrix.
 * @details Used e.g. for 0.5 * (K + trans(K)) symmetrization of a system
 * matrix. Unlike the dense overload this MATERIALIZES the result (one
 * O(nnz) copy): Eigen's lazy sparse transpose flips the storage order, and
 * sparse binary operations require both sides to share the same order.
 * @tparam TDerived Concrete type of the sparse matrix expression.
 * @param rM Sparse matrix to transpose.
 * @return The transposed matrix, materialized with the original storage order.
 */
template<class TDerived>
inline typename TDerived::PlainObject trans(const Eigen::SparseMatrixBase<TDerived>& rM)
{
    return rM.transpose();
}

/**
 * @brief 1-norm (sum of absolute values).
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rX Vector or matrix to compute the norm of.
 * @return The 1-norm of rX.
 */
template<class TDerived>
inline auto norm_1(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.template lpNorm<1>();
}

/**
 * @brief Euclidean norm (Frobenius norm for matrices, as in ublas).
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rX Vector or matrix to compute the norm of.
 * @return The Euclidean/Frobenius norm of rX.
 */
template<class TDerived>
inline auto norm_2(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.norm();
}

/**
 * @brief Squared Euclidean norm (avoids the sqrt/rsqrt of norm_2, as in ublas).
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rX Vector or matrix to compute the squared norm of.
 * @return The squared Euclidean norm of rX.
 */
template<class TDerived>
inline auto norm_2_square(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.squaredNorm();
}

/**
 * @brief Infinity norm (maximum absolute value).
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rX Vector or matrix to compute the norm of.
 * @return The infinity norm of rX.
 */
template<class TDerived>
inline auto norm_inf(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.template lpNorm<Eigen::Infinity>();
}

/**
 * @brief Sum of all coefficients.
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rX Vector or matrix whose coefficients are summed.
 * @return The sum of all coefficients of rX.
 */
template<class TDerived>
inline auto sum(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.sum();
}

/**
 * @brief Frobenius norm (same as norm_2 for vectors, as in ublas).
 * @tparam TDerived Concrete type of the Eigen expression.
 * @param rM Vector or matrix to compute the norm of.
 * @return The Frobenius norm of rM.
 */
template<class TDerived>
inline auto norm_frobenius(const Eigen::MatrixBase<TDerived>& rM)
{
    return rM.norm();
}

/**
 * @brief Element-wise (Hadamard) product.
 * @tparam TDerived1 Concrete type of the first operand.
 * @tparam TDerived2 Concrete type of the second operand.
 * @param rX First operand.
 * @param rY Second operand.
 * @return The element-wise product as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
inline auto element_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.cwiseProduct(rY.derived());
}

/**
 * @brief Element-wise division.
 * @tparam TDerived1 Concrete type of the numerator.
 * @tparam TDerived2 Concrete type of the denominator.
 * @param rX Numerator.
 * @param rY Denominator.
 * @return The element-wise quotient as a lazy Eigen expression.
 */
template<class TDerived1, class TDerived2>
inline auto element_div(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.cwiseQuotient(rY.derived());
}

namespace Internals {

/**
 * @class EigenDynamicNoAliasProxy
 * @brief noalias proxy for the dynamic Eigen-backed types.
 * @details Eigen RHS go through Eigen's own NoAlias (the exact ublas noalias
 * contract: assignment without a protective temporary); uBLAS expression RHS
 * (ZeroMatrix, IdentityMatrix, prod(...) of uBLAS operands, ...) are
 * evaluated element-wise, resizing the target as plain ublas assignment would.
 * @tparam TTarget Concrete Eigen-backed target type (EigenMatrix<T> or EigenVector<T>).
 */
template<class TTarget>
class EigenDynamicNoAliasProxy
{
public:

    /**
     * @brief Constructor wrapping the target that operator=/+=/-= write into.
     * @param rTarget Target to be assigned into without a protective temporary.
     */
    explicit EigenDynamicNoAliasProxy(TTarget& rTarget) : mrTarget(rTarget) {}

    /**
     * @brief Assignment from an Eigen expression, through Eigen's own NoAlias.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rExpression Source Eigen expression.
     * @return Reference to the wrapped target.
     */
    template<class TDerived>
    TTarget& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() = rExpression.derived();
        return mrTarget;
    }

    /**
     * @brief In-place addition of an Eigen expression, through Eigen's own NoAlias.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rExpression Source Eigen expression.
     * @return Reference to the wrapped target.
     */
    template<class TDerived>
    TTarget& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() += rExpression.derived();
        return mrTarget;
    }

    /**
     * @brief In-place subtraction of an Eigen expression, through Eigen's own NoAlias.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rExpression Source Eigen expression.
     * @return Reference to the wrapped target.
     */
    template<class TDerived>
    TTarget& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() -= rExpression.derived();
        return mrTarget;
    }

    /**
     * @brief Assignment from a dense uBLAS matrix expression, evaluated element-wise.
     * @details Resizes the target to match rExpression first, as plain ublas
     * assignment would.
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     * @return Reference to the wrapped target.
     */
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

    /**
     * @brief In-place addition of a dense uBLAS matrix expression, evaluated element-wise.
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     * @return Reference to the wrapped target.
     */
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                mrTarget(i, j) += r_e(i, j);
        return mrTarget;
    }

    /**
     * @brief In-place subtraction of a dense uBLAS matrix expression, evaluated element-wise.
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     * @return Reference to the wrapped target.
     */
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                mrTarget(i, j) -= r_e(i, j);
        return mrTarget;
    }

    /**
     * @brief Assignment from a dense uBLAS vector expression, evaluated element-wise.
     * @details Resizes the target to match rExpression first, as plain ublas
     * assignment would.
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     * @return Reference to the wrapped target.
     */
    template<class TExpression>
    TTarget& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        mrTarget.resize(r_e.size(), false);
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] = r_e(i);
        return mrTarget;
    }

    /**
     * @brief In-place addition of a dense uBLAS vector expression, evaluated element-wise.
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     * @return Reference to the wrapped target.
     */
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] += r_e(i);
        return mrTarget;
    }

    /**
     * @brief In-place subtraction of a dense uBLAS vector expression, evaluated element-wise.
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     * @return Reference to the wrapped target.
     */
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) mrTarget[i] -= r_e(i);
        return mrTarget;
    }

private:
    TTarget& mrTarget; /// The wrapped target, written into by operator=/+=/-=.
};

} // namespace Internals

/**
 * @brief noalias proxy accepting both Eigen and uBLAS right-hand sides, for an EigenMatrix target.
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Target matrix to be assigned into without a protective temporary.
 * @return A proxy through which operator=/+=/-= write into rM.
 */
template<class TDataType>
inline auto noalias(EigenMatrix<TDataType>& rM)
{
    return Internals::EigenDynamicNoAliasProxy<EigenMatrix<TDataType>>(rM);
}

/**
 * @brief noalias proxy accepting both Eigen and uBLAS right-hand sides, for an EigenVector target.
 * @tparam TDataType Scalar type of the vector.
 * @param rV Target vector to be assigned into without a protective temporary.
 * @return A proxy through which operator=/+=/-= write into rV.
 */
template<class TDataType>
inline auto noalias(EigenVector<TDataType>& rV)
{
    return Internals::EigenDynamicNoAliasProxy<EigenVector<TDataType>>(rV);
}

/**
 * @brief Row proxy (readable and writable, as ublas row()).
 * @details The uBLAS matrix_row models a *vector*; Eigen's vector convention
 * is a column, so the row block is transposed to keep vector semantics
 * (assignment to array_1d/EigenVector, inner_prod, v -= row(M,i), ...).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the row from.
 * @param I Zero-based row index.
 * @return The I-th row of rM, transposed into a column-shaped, writable Eigen expression.
 */
template<class TDataType>
inline auto row(EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
}

/**
 * @brief Row proxy (const version, as ublas row()).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the row from.
 * @param I Zero-based row index.
 * @return The I-th row of rM, transposed into a column-shaped, read-only Eigen expression.
 */
template<class TDataType>
inline auto row(const EigenMatrix<TDataType>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
}

/**
 * @brief Column proxy (readable and writable, as ublas column()).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the column from.
 * @param J Zero-based column index.
 * @return The J-th column of rM as a writable Eigen expression.
 */
template<class TDataType>
inline auto column(EigenMatrix<TDataType>& rM, const std::size_t J)
{
    return rM.col(J);
}

/**
 * @brief Column proxy (const version, as ublas column()).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the column from.
 * @param J Zero-based column index.
 * @return The J-th column of rM as a read-only Eigen expression.
 */
template<class TDataType>
inline auto column(const EigenMatrix<TDataType>& rM, const std::size_t J)
{
    return rM.col(J);
}

/**
 * @brief Vector subrange [Low, High) proxy (readable and writable, as ublas subrange()).
 * @tparam TDataType Scalar type of the vector.
 * @param rV Vector to extract the subrange from.
 * @param Low Zero-based index of the first entry included in the subrange.
 * @param High Zero-based index one past the last entry included in the subrange.
 * @return The [Low, High) subrange of rV as a writable Eigen expression.
 */
template<class TDataType>
inline auto subrange(EigenVector<TDataType>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

/**
 * @brief Vector subrange [Low, High) proxy (const version, as ublas subrange()).
 * @tparam TDataType Scalar type of the vector.
 * @param rV Vector to extract the subrange from.
 * @param Low Zero-based index of the first entry included in the subrange.
 * @param High Zero-based index one past the last entry included in the subrange.
 * @return The [Low, High) subrange of rV as a read-only Eigen expression.
 */
template<class TDataType>
inline auto subrange(const EigenVector<TDataType>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

/**
 * @brief Matrix subrange [Row1, Row2) x [Col1, Col2) block (readable and writable, as ublas subrange()).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the block from.
 * @param Row1 Zero-based index of the first row included in the block.
 * @param Row2 Zero-based index one past the last row included in the block.
 * @param Col1 Zero-based index of the first column included in the block.
 * @param Col2 Zero-based index one past the last column included in the block.
 * @return The [Row1, Row2) x [Col1, Col2) block of rM as a writable Eigen expression.
 */
template<class TDataType>
inline auto subrange(EigenMatrix<TDataType>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

/**
 * @brief Matrix subrange [Row1, Row2) x [Col1, Col2) block (const version, as ublas subrange()).
 * @tparam TDataType Scalar type of the matrix.
 * @param rM Matrix to extract the block from.
 * @param Row1 Zero-based index of the first row included in the block.
 * @param Row2 Zero-based index one past the last row included in the block.
 * @param Col1 Zero-based index of the first column included in the block.
 * @param Col2 Zero-based index one past the last column included in the block.
 * @return The [Row1, Row2) x [Col1, Col2) block of rM as a read-only Eigen expression.
 */
template<class TDataType>
inline auto subrange(const EigenMatrix<TDataType>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

} // namespace Kratos
