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
//                   Riccardo Rossi
//

#pragma once

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "includes/eigen_dense_types.h"

/**
 * @brief uBLAS-style free functions and proxies for the Eigen types.
 * @details Everything generic Kratos code writes against the uBLAS idiom
 * (prod, inner_prod, outer_prod, trans, noalias, row, column, subrange,
 * project, norms, lu_factorize, ...) is provided here in namespace Kratos with
 * a pure-Eigen implementation and the uBLAS SEMANTICS (trans() of a vector is
 * the identity, prod(v, M) is M^T v, row() is a vector, ...). The functions
 * deduce on Eigen::MatrixBase / Eigen::SparseMatrixBase, so every Kratos dense
 * type (the dynamic and bounded wrappers, array_1d, the proxies below) and
 * every Eigen expression is accepted.
 */
namespace Kratos
{

namespace Internals
{

template<class TDerived> std::true_type IsEigenDenseImpl(const Eigen::MatrixBase<TDerived>&);
std::false_type IsEigenDenseImpl(...);
template<class TDerived> std::true_type IsEigenSparseImpl(const Eigen::SparseMatrixBase<TDerived>&);
std::false_type IsEigenSparseImpl(...);

/// True for every type deriving from Eigen::MatrixBase (dense Kratos types and Eigen expressions).
template<class T>
inline constexpr bool IsEigenDense = decltype(IsEigenDenseImpl(std::declval<const T&>()))::value;

/// True for every type deriving from Eigen::SparseMatrixBase.
template<class T>
inline constexpr bool IsEigenSparse = decltype(IsEigenSparseImpl(std::declval<const T&>()))::value;

/// True for the plain (storage-owning) Eigen objects.
template<class TDerived>
inline constexpr bool IsPlainObject = std::is_base_of_v<Eigen::PlainObjectBase<TDerived>, TDerived>;

} // namespace Internals

///@name Index ranges
///@{

/**
 * @class range
 * @brief Half-open index range [start, stop), as boost::numeric::ublas::range.
 */
class range
{
public:
    using size_type = std::size_t;

    range(const size_type Start, const size_type Stop) : mStart(Start), mSize(Stop - Start) {}

    size_type start() const { return mStart; }
    size_type size() const { return mSize; }
    bool empty() const { return mSize == 0; }
    size_type operator()(const size_type Index) const { return mStart + Index; }
    size_type operator[](const size_type Index) const { return mStart + Index; }

    /// Sentinel selecting the whole extent of the target (resolved by preprocess()).
    static range all() { return range(0, std::numeric_limits<size_type>::max()); }

    /// The range resolved against the extent of the target it is applied to.
    range preprocess(const size_type Extent) const
    {
        return (mStart == 0 && mSize == std::numeric_limits<size_type>::max()) ? range(0, Extent) : *this;
    }

    bool operator==(const range& rOther) const { return mStart == rOther.mStart && mSize == rOther.mSize; }

private:
    size_type mStart;
    size_type mSize;
};

/**
 * @class slice
 * @brief Strided index set start + stride * i, i in [0, size), as boost::numeric::ublas::slice.
 */
class slice
{
public:
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    slice(const size_type Start, const difference_type Stride, const size_type Size) : mStart(Start), mStride(Stride), mSize(Size) {}

    size_type start() const { return mStart; }
    difference_type stride() const { return mStride; }
    size_type size() const { return mSize; }
    bool empty() const { return mSize == 0; }
    size_type operator()(const size_type Index) const { return mStart + static_cast<size_type>(mStride * static_cast<difference_type>(Index)); }
    size_type operator[](const size_type Index) const { return (*this)(Index); }

    /// Sentinel selecting the whole extent of the target (resolved by preprocess()).
    static slice all() { return slice(0, 1, std::numeric_limits<size_type>::max()); }

    /// The slice resolved against the extent of the target it is applied to.
    slice preprocess(const size_type Extent) const
    {
        return (mStart == 0 && mStride == 1 && mSize == std::numeric_limits<size_type>::max()) ? slice(0, 1, Extent) : *this;
    }

private:
    size_type mStart;
    difference_type mStride;
    size_type mSize;
};

///@}
///@name Proxies
///@{
// The uBLAS proxies (matrix_row, matrix_column, vector_range, matrix_range,
// vector_slice) are thin subclasses of the corresponding Eigen block
// expressions: every Eigen operation applies to them (they ARE Eigen
// expressions, referencing the parent's storage), while the subclass adds the
// uBLAS constructors and the unsigned size()/size1()/size2() of the uBLAS
// surface. A matrix row is a VECTOR (column-shaped), as in uBLAS.

namespace Internals
{
template<class TExpression> using MatrixRowBase = std::remove_cv_t<decltype(std::declval<TExpression&>().row(0).transpose())>;
template<class TExpression> using MatrixColumnBase = std::remove_cv_t<decltype(std::declval<TExpression&>().col(0))>;
template<class TExpression> using VectorRangeBase = std::remove_cv_t<decltype(std::declval<TExpression&>().segment(0, 0))>;
template<class TExpression> using MatrixRangeBase = std::remove_cv_t<decltype(std::declval<TExpression&>().block(0, 0, 0, 0))>;
template<class TExpression> using VectorSliceBase = Eigen::Map<
    std::conditional_t<std::is_const_v<TExpression>, const DynamicColumnVector<typename TExpression::value_type>, DynamicColumnVector<typename TExpression::value_type>>,
    Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>;
} // namespace Internals

/// Row proxy (boost::numeric::ublas::matrix_row): a column-shaped, writable view of one matrix row.
template<class TExpression>
class matrix_row : public Internals::MatrixRowBase<TExpression>
{
public:
    using BaseType = Internals::MatrixRowBase<TExpression>;
    using value_type = typename BaseType::Scalar;
    using size_type = std::size_t;

    matrix_row(TExpression& rMatrix, const std::size_t Index) : BaseType(rMatrix.row(static_cast<Eigen::Index>(Index)).transpose()) {}
    matrix_row(const BaseType& rExpression) : BaseType(rExpression) {}
    matrix_row(const matrix_row& rOther) = default;

    matrix_row& operator=(const matrix_row& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /// uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view.
    template<class TDerived> matrix_row& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_row& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_row& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/// Column proxy (boost::numeric::ublas::matrix_column): a writable view of one matrix column.
template<class TExpression>
class matrix_column : public Internals::MatrixColumnBase<TExpression>
{
public:
    using BaseType = Internals::MatrixColumnBase<TExpression>;
    using value_type = typename BaseType::Scalar;
    using size_type = std::size_t;

    matrix_column(TExpression& rMatrix, const std::size_t Index) : BaseType(rMatrix.col(static_cast<Eigen::Index>(Index))) {}
    matrix_column(const BaseType& rExpression) : BaseType(rExpression) {}
    matrix_column(const matrix_column& rOther) = default;

    matrix_column& operator=(const matrix_column& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /// uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view.
    template<class TDerived> matrix_column& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_column& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_column& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/// Range proxy (boost::numeric::ublas::vector_range): a writable view of a contiguous vector subrange.
template<class TExpression>
class vector_range : public Internals::VectorRangeBase<TExpression>
{
public:
    using BaseType = Internals::VectorRangeBase<TExpression>;
    using value_type = typename BaseType::Scalar;
    using size_type = std::size_t;

    vector_range(TExpression& rVector, const range& rRange)
        : BaseType(rVector.segment(static_cast<Eigen::Index>(rRange.preprocess(rVector.size()).start()), static_cast<Eigen::Index>(rRange.preprocess(rVector.size()).size()))) {}
    vector_range(TExpression& rVector, const std::size_t Start, const std::size_t Stop) : vector_range(rVector, range(Start, Stop)) {}
    vector_range(const BaseType& rExpression) : BaseType(rExpression) {}
    vector_range(const vector_range& rOther) = default;

    vector_range& operator=(const vector_range& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /// uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view.
    template<class TDerived> vector_range& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> vector_range& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> vector_range& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/// Range proxy (boost::numeric::ublas::matrix_range): a writable view of a matrix block.
template<class TExpression>
class matrix_range : public Internals::MatrixRangeBase<TExpression>
{
public:
    using BaseType = Internals::MatrixRangeBase<TExpression>;
    using value_type = typename BaseType::Scalar;
    using size_type = std::size_t;

    matrix_range(TExpression& rMatrix, const range& rRange1, const range& rRange2)
        : BaseType(rMatrix.block(
            static_cast<Eigen::Index>(rRange1.preprocess(rMatrix.rows()).start()),
            static_cast<Eigen::Index>(rRange2.preprocess(rMatrix.cols()).start()),
            static_cast<Eigen::Index>(rRange1.preprocess(rMatrix.rows()).size()),
            static_cast<Eigen::Index>(rRange2.preprocess(rMatrix.cols()).size()))) {}
    matrix_range(const BaseType& rExpression) : BaseType(rExpression) {}
    matrix_range(const matrix_range& rOther) = default;

    matrix_range& operator=(const matrix_range& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /// uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view.
    template<class TDerived> matrix_range& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_range& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> matrix_range& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    std::size_t size1() const { return static_cast<std::size_t>(BaseType::rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(BaseType::cols()); }
};

/// Slice proxy (boost::numeric::ublas::vector_slice): a writable strided view
/// over the contiguous storage of a vector.
template<class TExpression>
class vector_slice : public Internals::VectorSliceBase<TExpression>
{
public:
    using BaseType = Internals::VectorSliceBase<TExpression>;
    using value_type = typename BaseType::Scalar;
    using size_type = std::size_t;

    vector_slice(TExpression& rVector, const slice& rSlice)
        : BaseType(rVector.data().begin() + rSlice.preprocess(rVector.size()).start(),
                   static_cast<Eigen::Index>(rSlice.preprocess(rVector.size()).size()),
                   Eigen::InnerStride<Eigen::Dynamic>(static_cast<Eigen::Index>(rSlice.preprocess(rVector.size()).stride()))) {}
    vector_slice(const BaseType& rExpression) : BaseType(rExpression) {}
    vector_slice(const vector_slice& rOther) = default;

    vector_slice& operator=(const vector_slice& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /// uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view.
    template<class TDerived> vector_slice& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> vector_slice& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }
    template<class TDerived> vector_slice& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/// Row of a matrix as a (column-shaped) vector proxy, readable and writable (ublas row()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_row<TMatrix> row(TMatrix& rM, const std::size_t I) { return matrix_row<TMatrix>(rM, I); }

/// Row of a matrix as a (column-shaped) vector proxy, read-only (ublas row()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_row<const TMatrix> row(const TMatrix& rM, const std::size_t I) { return matrix_row<const TMatrix>(rM, I); }

/// Column of a matrix as a vector proxy, readable and writable (ublas column()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_column<TMatrix> column(TMatrix& rM, const std::size_t J) { return matrix_column<TMatrix>(rM, J); }

/// Column of a matrix as a vector proxy, read-only (ublas column()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_column<const TMatrix> column(const TMatrix& rM, const std::size_t J) { return matrix_column<const TMatrix>(rM, J); }

/// Vector subrange [Low, High), readable and writable (ublas subrange()).
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<TVector> subrange(TVector& rV, const std::size_t Low, const std::size_t High) { return vector_range<TVector>(rV, Low, High); }

/// Vector subrange [Low, High), read-only (ublas subrange()).
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<const TVector> subrange(const TVector& rV, const std::size_t Low, const std::size_t High) { return vector_range<const TVector>(rV, Low, High); }

/// Matrix block [Row1, Row2) x [Col1, Col2), readable and writable (ublas subrange()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<TMatrix> subrange(TMatrix& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return matrix_range<TMatrix>(rM, range(Row1, Row2), range(Col1, Col2));
}

/// Matrix block [Row1, Row2) x [Col1, Col2), read-only (ublas subrange()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<const TMatrix> subrange(const TMatrix& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return matrix_range<const TMatrix>(rM, range(Row1, Row2), range(Col1, Col2));
}

/// Vector range projection, readable and writable (ublas project()).
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<TVector> project(TVector& rV, const range& rRange) { return vector_range<TVector>(rV, rRange); }

/// Vector range projection, read-only (ublas project()).
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<const TVector> project(const TVector& rV, const range& rRange) { return vector_range<const TVector>(rV, rRange); }

/// Matrix range projection, readable and writable (ublas project()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<TMatrix> project(TMatrix& rM, const range& rRange1, const range& rRange2) { return matrix_range<TMatrix>(rM, rRange1, rRange2); }

/// Matrix range projection, read-only (ublas project()).
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<const TMatrix> project(const TMatrix& rM, const range& rRange1, const range& rRange2) { return matrix_range<const TMatrix>(rM, rRange1, rRange2); }

/// Vector slice projection over a storage-owning vector, readable and writable (ublas project()).
template<class TVector> requires (Internals::IsEigenDense<TVector> && Internals::IsPlainObject<typename TVector::BaseType>)
inline vector_slice<TVector> project(TVector& rV, const slice& rSlice) { return vector_slice<TVector>(rV, rSlice); }

/// Vector slice projection over a storage-owning vector, read-only (ublas project()).
template<class TVector> requires (Internals::IsEigenDense<TVector> && Internals::IsPlainObject<typename TVector::BaseType>)
inline vector_slice<const TVector> project(const TVector& rV, const slice& rSlice) { return vector_slice<const TVector>(rV, rSlice); }

///@}
///@name Products
///@{

namespace Internals
{
/**
 * @brief Shared implementation for the dense prod() overloads below.
 * @details Kept out of the "prod" overload set on purpose: a call from inside
 * a prod() overload to unqualified prod(rA, rB) would re-run overload
 * resolution against every prod overload (including the sparse ones), and
 * deducing those against a reference-to-MatrixBase argument forces an
 * instantiation of Eigen::SparseMatrixBase<TDerived> that hard-errors.
 * A vector first operand keeps the ublas prod(v, M) semantics (v^T M),
 * returned column-shaped so it assigns to the vector types.
 */
template<class TDerived1, class TDerived2>
inline auto dense_prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    if constexpr (IsColumnShaped<TDerived1> && TDerived2::RowsAtCompileTime != 1) {
        // prod(v, M) = v^T M: an N x 1 first operand is a vector by type in
        // uBLAS unless the second operand is a 1 x K matrix (then the N x 1
        // operand is a bounded matrix and this is the plain N x K product).
        return (rA.derived().transpose() * rB.derived()).transpose();
    } else {
        return rA.derived() * rB.derived();
    }
}
} // namespace Internals

/**
 * @brief Matrix/matrix and matrix/vector product (lazy Eigen expression).
 * @details The requires-clause guards TDerived1/TDerived2 the same way as the
 * sparse overloads below: prod<TResult>(...) shares the "prod" name, so a call
 * like prod<Matrix>(...) explicitly substitutes TDerived1 = Matrix into every
 * prod overload, and naming Eigen::MatrixBase<Matrix> for that substitution
 * forces Eigen to complete it, which hard-errors (Eigen::internal::traits is
 * unspecialized for the Kratos wrapper types). The leading TShield parameter
 * (never specified by callers) absorbs that explicit argument.
 */
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires {
        typename Eigen::internal::traits<TDerived1>::StorageKind;
        typename Eigen::internal::traits<TDerived2>::StorageKind;
    })
inline auto prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return Internals::dense_prod(rA, rB);
}

/// Product with an explicitly pinned result type, as ublas prod<Vector>(A, B).
template<class TResult, class TDerived1, class TDerived2>
inline TResult prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return TResult(Internals::dense_prod(rA, rB));
}

/// Sparse-matrix/vector (or sparse/dense) product.
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires { typename Eigen::internal::traits<TDerived1>::StorageKind; })
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rX)
{
    return rA.derived() * rX.derived();
}

/// Vector/sparse-matrix product: v^T A, returned column-shaped (uBLAS prod(v, A)).
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires { typename Eigen::internal::traits<TDerived2>::StorageKind; })
inline auto prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::SparseMatrixBase<TDerived2>& rA)
{
    return rA.derived().transpose() * rX.derived();
}

/// Sparse-matrix/sparse-matrix product.
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires {
        typename Eigen::internal::traits<TDerived1>::StorageKind;
        typename Eigen::internal::traits<TDerived2>::StorageKind;
    })
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::SparseMatrixBase<TDerived2>& rB)
{
    return rA.derived() * rB.derived();
}

/// Scalar (dot) product of two vectors. As in uBLAS the product is NOT
/// conjugated for complex scalars (Eigen's dot() would conjugate).
template<class TDerived1, class TDerived2>
inline auto inner_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    using Scalar1 = typename TDerived1::Scalar;
    using Scalar2 = typename TDerived2::Scalar;
    if constexpr (std::is_same_v<Scalar1, Scalar2>) {
        return (rX.derived().transpose() * rY.derived()).value();
    } else {
        using ResultScalar = decltype(std::declval<Scalar1>() * std::declval<Scalar2>());
        return (rX.derived().template cast<ResultScalar>().transpose() * rY.derived().template cast<ResultScalar>()).value();
    }
}

/// Outer product of two vectors, rX * rY^T; a row-shaped rY (as produced by
/// Eigen's own transpose) is accepted as is.
template<class TDerived1, class TDerived2>
inline auto outer_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    if constexpr (TDerived2::RowsAtCompileTime == 1 && TDerived2::ColsAtCompileTime != 1) {
        return rX.derived() * rY.derived();
    } else {
        return rX.derived() * rY.derived().transpose();
    }
}

/**
 * @brief uBLAS-style axpy product: rY = A B (Init = true) or rY += A B.
 * @details A single overload covers the three uBLAS forms, dispatching on the
 * compile-time shape of the first operand exactly as prod() does: matrix x
 * vector, matrix x matrix and, for a vector first operand, the vector x
 * matrix form (y = M^T v). As in uBLAS the target must not alias either
 * operand; the assignment goes through Eigen's noalias() path.
 */
template<class TDerived1, class TDerived2, class TTargetType>
inline TTargetType& axpy_prod(
    const Eigen::MatrixBase<TDerived1>& rA,
    const Eigen::MatrixBase<TDerived2>& rB,
    TTargetType& rY,
    const bool Init = true)
{
    using Scalar = typename TTargetType::Scalar;
    if constexpr (Internals::IsColumnShaped<TDerived1> && TDerived2::RowsAtCompileTime != 1) {
        if (Init) rY.noalias() = (rB.derived().transpose() * rA.derived()).template cast<Scalar>();
        else rY.noalias() += (rB.derived().transpose() * rA.derived()).template cast<Scalar>();
    } else {
        if (Init) rY.noalias() = (rA.derived() * rB.derived()).template cast<Scalar>();
        else rY.noalias() += (rA.derived() * rB.derived()).template cast<Scalar>();
    }
    return rY;
}

/// Sparse counterpart of the axpy product (sparse matrix times a dense vector or matrix).
template<class TDerived1, class TDerived2, class TTargetType>
inline TTargetType& axpy_prod(
    const Eigen::SparseMatrixBase<TDerived1>& rA,
    const Eigen::MatrixBase<TDerived2>& rX,
    TTargetType& rY,
    const bool Init = true)
{
    if (Init) rY.noalias() = rA.derived() * rX.derived();
    else rY.noalias() += rA.derived() * rX.derived();
    return rY;
}

///@}
///@name Transposition
///@{

/**
 * @brief Transpose, with the uBLAS TYPE-based semantics.
 * @details uBLAS defines trans() on a vector as the identity ((trans v)[i] = v[i],
 * boost vector_expression.hpp) and discriminates vector from matrix by type,
 * not by shape. Every compile-time vector-shaped expression (Vector, array_1d,
 * BoundedVector, the row/column proxies, ...) therefore comes back unchanged
 * (by reference for the storage-owning types, by value for lightweight
 * expressions), so idioms like outer_prod(x, trans(y)) and prod(trans(v), M)
 * mean the same under both backends; matrices are lazily transposed.
 */
template<class TDerived>
inline decltype(auto) trans(const Eigen::MatrixBase<TDerived>& rX)
{
    if constexpr (Internals::IsColumnShaped<TDerived>) {
        if constexpr (Internals::IsPlainObject<TDerived>) {
            return (rX.derived());
        } else {
            return rX.derived();
        }
    } else {
        return rX.transpose();
    }
}

/// A bounded matrix is a matrix by type whatever its static shape: an N x 1 or
/// 1 x N BoundedMatrix transposes (unlike the N-vectors, see above). This
/// overload is an exact match for the wrapper type, so it beats the generic
/// MatrixBase template.
template<class T, std::size_t TSize1, std::size_t TSize2>
inline auto trans(const EigenBoundedMatrix<T, TSize1, TSize2>& rM)
{
    return rM.transpose();
}

/// Transpose of a sparse matrix. Unlike the dense overload this MATERIALIZES
/// the result: Eigen's lazy sparse transpose flips the storage order, and
/// sparse binary operations require both sides to share the same order.
template<class TDerived>
inline typename TDerived::PlainObject trans(const Eigen::SparseMatrixBase<TDerived>& rM)
{
    return rM.transpose();
}

///@}
///@name Norms and reductions
///@{

/// 1-norm: the sum of absolute values for a vector, the maximum absolute
/// column sum for a matrix (ublas semantics; 0 for an empty matrix).
template<class TDerived>
inline typename Eigen::NumTraits<typename TDerived::Scalar>::Real norm_1(const Eigen::MatrixBase<TDerived>& rX)
{
    if constexpr (Internals::IsColumnShaped<TDerived>) {
        return rX.template lpNorm<1>();
    } else {
        if (rX.size() == 0) return 0; // maxCoeff() asserts on an empty expression, ublas returns 0
        return rX.cwiseAbs().colwise().sum().maxCoeff();
    }
}

/// Euclidean norm (Frobenius norm for matrices, as in ublas).
template<class TDerived>
inline auto norm_2(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.norm();
}

/// Squared Euclidean norm.
template<class TDerived>
inline auto norm_2_square(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.squaredNorm();
}

/// Infinity norm: the maximum absolute value for a vector, the maximum
/// absolute row sum for a matrix (ublas semantics; 0 for an empty matrix).
template<class TDerived>
inline typename Eigen::NumTraits<typename TDerived::Scalar>::Real norm_inf(const Eigen::MatrixBase<TDerived>& rX)
{
    if constexpr (Internals::IsColumnShaped<TDerived>) {
        return rX.template lpNorm<Eigen::Infinity>();
    } else {
        if (rX.size() == 0) return 0; // maxCoeff() asserts on an empty expression, ublas returns 0
        return rX.cwiseAbs().rowwise().sum().maxCoeff();
    }
}

/// Index of the first entry of maximum absolute value (ublas index_norm_inf).
template<class TDerived>
inline std::size_t index_norm_inf(const Eigen::MatrixBase<TDerived>& rX)
{
    using std::abs;
    std::size_t index = 0;
    auto max_abs = abs(rX.derived().coeff(0));
    for (Eigen::Index i = 1; i < rX.size(); ++i) {
        const auto a = abs(rX.derived().coeff(i));
        if (a > max_abs) {
            max_abs = a;
            index = static_cast<std::size_t>(i);
        }
    }
    return index;
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

/// Element-wise (Hadamard) product.
template<class TDerived1, class TDerived2>
inline auto element_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.derived().cwiseProduct(rY.derived());
}

/// Element-wise division.
template<class TDerived1, class TDerived2>
inline auto element_div(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.derived().cwiseQuotient(rY.derived());
}

///@}
///@name noalias
///@{

namespace Internals
{
/**
 * @class NoAliasProxy
 * @brief noalias() proxy: assignment without a protective temporary, through
 * Eigen's own NoAlias, with the ublas element-wise scalar conversion (the
 * cast is the identity for matching scalars).
 * @tparam TTargetHolder The target expression type, held by reference (lvalue
 * targets) or by value (temporary proxies such as row(M, i)).
 */
template<class TTargetHolder>
class NoAliasProxy
{
public:
    using TargetType = std::remove_reference_t<TTargetHolder>;
    using Scalar = typename TargetType::Scalar;

    explicit NoAliasProxy(TargetType& rTarget) : mTarget(rTarget) {}
    explicit NoAliasProxy(TargetType&& rTarget) : mTarget(std::move(rTarget)) {}

    template<class TDerived>
    TargetType& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() = rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    template<class TDerived>
    TargetType& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() += rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    template<class TDerived>
    TargetType& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() -= rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    template<class TDerived>
    TargetType& operator=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget = rExpression.derived();
        return mTarget;
    }

    template<class TDerived>
    TargetType& operator+=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget += rExpression.derived();
        return mTarget;
    }

    template<class TDerived>
    TargetType& operator-=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget -= rExpression.derived();
        return mTarget;
    }

private:
    TTargetHolder mTarget;
};
} // namespace Internals

/// noalias() on a dense lvalue target (a container or a stored proxy).
template<class TDerived>
inline auto noalias(Eigen::MatrixBase<TDerived>& rTarget)
{
    return Internals::NoAliasProxy<TDerived&>(rTarget.derived());
}

/// noalias() on a temporary dense proxy, e.g. noalias(row(M, i)) = ...
template<class TDerived>
inline auto noalias(Eigen::MatrixBase<TDerived>&& rTarget)
{
    return Internals::NoAliasProxy<TDerived>(std::move(rTarget.derived()));
}

/// noalias() on a sparse lvalue target.
template<class TDerived>
inline auto noalias(Eigen::SparseMatrixBase<TDerived>& rTarget)
{
    return Internals::NoAliasProxy<TDerived&>(rTarget.derived());
}

///@}
///@name LU factorization (boost::numeric::ublas::lu_factorize / lu_substitute)
///@{

/// Permutation (transposition sequence) produced by lu_factorize, as boost::numeric::ublas::permutation_matrix.
template<class T = std::size_t>
class permutation_matrix : public EigenVector<T>
{
public:
    using BaseType = EigenVector<T>;
    using size_type = std::size_t;

    explicit permutation_matrix(const std::size_t Size) : BaseType(Size)
    {
        for (std::size_t i = 0; i < Size; ++i) {
            (*this)[i] = static_cast<T>(i);
        }
    }

    /// Applies the transposition sequence to the rows of a vector or matrix (ublas swap_rows).
    template<class TDerived>
    void ApplyTo(Eigen::MatrixBase<TDerived>& rX) const
    {
        for (std::size_t i = 0; i < this->size(); ++i) {
            const std::size_t pivot = static_cast<std::size_t>((*this)[i]);
            if (pivot != i) {
                rX.row(static_cast<Eigen::Index>(i)).swap(rX.row(static_cast<Eigen::Index>(pivot)));
            }
        }
    }
};

/**
 * @brief In-place LU factorization with partial pivoting, as boost::numeric::ublas::lu_factorize(m, pm).
 * @details A direct port of the uBLAS algorithm, so the packed L\U layout, the
 * transposition sequence written into rPM (used e.g. by the determinant sign
 * in MathUtils) and the return value (0, or i + 1 for the first zero pivot)
 * are identical to the uBLAS ones.
 */
template<class TDerived, class TPermutation>
inline std::size_t lu_factorize(Eigen::MatrixBase<TDerived>& rM, TPermutation& rPM)
{
    using Scalar = typename TDerived::Scalar;
    using std::abs;
    auto& r_m = rM.derived();
    const Eigen::Index size1 = r_m.rows();
    const Eigen::Index size2 = r_m.cols();
    const Eigen::Index size = std::min(size1, size2);
    std::size_t singular = 0;
    for (Eigen::Index i = 0; i < size; ++i) {
        // pivot: first row of maximum absolute value in column i, from the diagonal down
        Eigen::Index i_norm_inf = i;
        auto max_abs = abs(r_m(i, i));
        for (Eigen::Index k = i + 1; k < size1; ++k) {
            const auto a = abs(r_m(k, i));
            if (a > max_abs) {
                max_abs = a;
                i_norm_inf = k;
            }
        }
        if (r_m(i_norm_inf, i) != Scalar(0)) {
            if (i_norm_inf != i) {
                rPM(static_cast<std::size_t>(i)) = static_cast<std::size_t>(i_norm_inf);
                r_m.row(i).swap(r_m.row(i_norm_inf));
            }
            r_m.col(i).tail(size1 - i - 1) *= Scalar(1) / r_m(i, i);
        } else if (singular == 0) {
            singular = static_cast<std::size_t>(i + 1);
        }
        r_m.bottomRightCorner(size1 - i - 1, size2 - i - 1).noalias() -= r_m.col(i).tail(size1 - i - 1) * r_m.row(i).tail(size2 - i - 1);
    }
    return singular;
}

/// Forward/back substitution with the factors of lu_factorize, as
/// boost::numeric::ublas::lu_substitute(m, pm, x); rX (a vector or a matrix
/// of right-hand sides) is overwritten with the solution.
template<class TDerived, class TPermutation, class TDerivedRhs>
inline void lu_substitute(const Eigen::MatrixBase<TDerived>& rM, const TPermutation& rPM, Eigen::MatrixBase<TDerivedRhs>& rX)
{
    auto& r_x = rX.derived();
    for (std::size_t i = 0; i < static_cast<std::size_t>(rPM.size()); ++i) {
        const std::size_t pivot = static_cast<std::size_t>(rPM(i));
        if (pivot != i) {
            r_x.row(static_cast<Eigen::Index>(i)).swap(r_x.row(static_cast<Eigen::Index>(pivot)));
        }
    }
    rM.derived().template triangularView<Eigen::UnitLower>().solveInPlace(r_x);
    rM.derived().template triangularView<Eigen::Upper>().solveInPlace(r_x);
}

///@}
///@name Indirect access (boost::numeric::ublas::indirect_array / matrix_indirect)
///@{

/// Index array of an indirect view: any indexable container works (ublas indirect_array<A>).
template<class TArray>
using indirect_array = TArray;

namespace Internals
{
/// Nullary functor reading a matrix through two index arrays.
template<class TMatrix, class TIndexArray>
struct IndirectAccessOp
{
    const TMatrix* mpMatrix;
    const TIndexArray* mpRowIndices;
    const TIndexArray* mpColumnIndices;

    typename TMatrix::value_type operator()(const Eigen::Index I, const Eigen::Index J) const
    {
        return (*mpMatrix)((*mpRowIndices)[static_cast<std::size_t>(I)], (*mpColumnIndices)[static_cast<std::size_t>(J)]);
    }
};
} // namespace Internals

/// Read-only indirect (gather) view of a matrix through two index arrays, as
/// boost::numeric::ublas::matrix_indirect; a genuine Eigen expression, so it
/// can be evaluated into any dense type.
template<class TMatrix, class TIndexArray>
class matrix_indirect : public Eigen::CwiseNullaryOp<Internals::IndirectAccessOp<std::remove_const_t<TMatrix>, TIndexArray>, Internals::DynamicRowMajorMatrix<typename std::remove_const_t<TMatrix>::value_type>>
{
public:
    using OpType = Internals::IndirectAccessOp<std::remove_const_t<TMatrix>, TIndexArray>;
    using BaseType = Eigen::CwiseNullaryOp<OpType, Internals::DynamicRowMajorMatrix<typename std::remove_const_t<TMatrix>::value_type>>;
    using value_type = typename std::remove_const_t<TMatrix>::value_type;
    using size_type = std::size_t;

    matrix_indirect(const TMatrix& rMatrix, const TIndexArray& rRowIndices, const TIndexArray& rColumnIndices)
        : BaseType(static_cast<Eigen::Index>(rRowIndices.size()), static_cast<Eigen::Index>(rColumnIndices.size()), OpType{&rMatrix, &rRowIndices, &rColumnIndices}) {}

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

///@}

} // namespace Kratos
