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
 * @brief uBLAS-style free functions and proxies for the Eigen types
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

/**
 * @brief Detection overload selected for the types deriving from Eigen::MatrixBase
 * @details Declared but never defined, it is only used in an unevaluated context by IsEigenDense. Overload resolution
 * picks it when the argument derives from Eigen::MatrixBase<TDerived>, and it yields std::true_type.
 * @tparam TDerived Derived type of the Eigen::MatrixBase argument, deduced from the call.
 * @return std::true_type, the detection result when the argument derives from Eigen::MatrixBase.
 */
template<class TDerived> std::true_type IsEigenDenseImpl(const Eigen::MatrixBase<TDerived>&);

/**
 * @brief Fallback detection overload for the types not deriving from Eigen::MatrixBase
 * @details Declared but never defined, it is only used in an unevaluated context by IsEigenDense. The ellipsis makes it
 * the worst match in overload resolution, so it is selected only when the template overload is not viable, and it
 * yields std::false_type.
 * @return std::false_type, the detection result when the argument does not derive from Eigen::MatrixBase.
 */
std::false_type IsEigenDenseImpl(...);

/**
 * @brief Detection overload selected for the types deriving from Eigen::SparseMatrixBase
 * @details Declared but never defined, it is only used in an unevaluated context by IsEigenSparse. Overload resolution
 * picks it when the argument derives from Eigen::SparseMatrixBase<TDerived>, and it yields std::true_type.
 * @tparam TDerived Derived type of the Eigen::SparseMatrixBase argument, deduced from the call.
 * @return std::true_type, the detection result when the argument derives from Eigen::SparseMatrixBase.
 */
template<class TDerived> std::true_type IsEigenSparseImpl(const Eigen::SparseMatrixBase<TDerived>&);

/**
 * @brief Fallback detection overload for the types not deriving from Eigen::SparseMatrixBase
 * @details Declared but never defined, it is only used in an unevaluated context by IsEigenSparse. The ellipsis makes
 * it the worst match in overload resolution, so it is selected only when the template overload is not viable, and it
 * yields std::false_type.
 * @return std::false_type, the detection result when the argument does not derive from Eigen::SparseMatrixBase.
 */
std::false_type IsEigenSparseImpl(...);

/**
 * @brief True for every type deriving from Eigen::MatrixBase
 * @details This covers the dense Kratos types and the Eigen expressions. It is evaluated by overload resolution of
 * IsEigenDenseImpl on a const reference to T.
 * @tparam T Type tested for deriving from Eigen::MatrixBase.
 */
template<class T>
inline constexpr bool IsEigenDense = decltype(IsEigenDenseImpl(std::declval<const T&>()))::value;

/**
 * @brief True for every type deriving from Eigen::SparseMatrixBase
 * @details It is evaluated by overload resolution of IsEigenSparseImpl on a const reference to T.
 * @tparam T Type tested for deriving from Eigen::SparseMatrixBase.
 */
template<class T>
inline constexpr bool IsEigenSparse = decltype(IsEigenSparseImpl(std::declval<const T&>()))::value;

/**
 * @brief True for the plain (storage-owning) Eigen objects
 * @details It holds when TDerived derives from Eigen::PlainObjectBase<TDerived>, that is for an Eigen object owning its
 * storage rather than an expression, block or map referencing storage owned elsewhere.
 * @tparam TDerived Eigen type tested, which is also the derived type given to Eigen::PlainObjectBase.
 */
template<class TDerived>
inline constexpr bool IsPlainObject = std::is_base_of_v<Eigen::PlainObjectBase<TDerived>, TDerived>;

} // namespace Internals

///@name Index ranges
///@{

/**
 * @brief Half-open index range [start, stop), as boost::numeric::ublas::range
 * @details It stores the first index and the number of indices (stop - start). The stop is expected not to be lower
 * than the start, since the size is computed with unsigned arithmetic. The sentinel returned by all() selects the whole
 * extent of the target, and it is resolved by preprocess().
 */
class range
{
public:
    /// Type of the indices and sizes of the range. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Constructs the range [Start, Stop)
     * @details The size is stored as Stop - Start, with no check that Stop is not lower than Start.
     * @param Start First index of the range.
     * @param Stop End of the range, the first index that does not belong to it.
     */
    range(const size_type Start, const size_type Stop) : mStart(Start), mSize(Stop - Start) {}

    /**
     * @brief First index of the range
     * @details Returns the start given at construction.
     * @return First index of the range.
     */
    size_type start() const { return mStart; }

    /**
     * @brief Number of indices in the range
     * @details Returns the stored size, that is Stop - Start.
     * @return Number of indices in the range.
     */
    size_type size() const { return mSize; }

    /**
     * @brief Checks whether the range contains no index
     * @details Returns true when the size is zero.
     * @return True if the range contains no index, false otherwise.
     */
    bool empty() const { return mSize == 0; }

    /**
     * @brief Index of the target addressed by the Index-th element of the range
     * @details Returns start + Index, with no bound check.
     * @param Index Position inside the range, counted from its first element.
     * @return Index of the target, that is start() + Index.
     */
    size_type operator()(const size_type Index) const { return mStart + Index; }

    /**
     * @brief Index of the target addressed by the Index-th element of the range
     * @details Same as operator(): returns start + Index, with no bound check.
     * @param Index Position inside the range, counted from its first element.
     * @return Index of the target, that is start() + Index.
     */
    size_type operator[](const size_type Index) const { return mStart + Index; }

    /**
     * @brief Sentinel selecting the whole extent of the target
     * @details It is the range starting at 0 with the maximum representable size, and it is resolved against the actual
     * extent by preprocess().
     * @return Sentinel range starting at 0 with the maximum representable size.
     */
    static range all() { return range(0, std::numeric_limits<size_type>::max()); }

    /**
     * @brief Range resolved against the extent of the target it is applied to
     * @details Returns range(0, Extent) when this range is the all() sentinel, otherwise this range unchanged.
     * @param Extent Actual extent of the target the range is applied to.
     * @return range(0, Extent) when this range is the all() sentinel, otherwise this range unchanged.
     */
    range preprocess(const size_type Extent) const
    {
        return (mStart == 0 && mSize == std::numeric_limits<size_type>::max()) ? range(0, Extent) : *this;
    }

    /**
     * @brief Equality comparison
     * @details Two ranges are equal when they have the same start and the same size.
     * @param rOther Range to compare with.
     * @return True if both ranges have the same start and the same size, false otherwise.
     */
    bool operator==(const range& rOther) const { return mStart == rOther.mStart && mSize == rOther.mSize; }

private:
    /**
     * @brief First index of the range
     * @details Set at construction.
     */
    size_type mStart;

    /**
     * @brief Number of indices in the range
     * @details Computed at construction as Stop - Start.
     */
    size_type mSize;
};

/**
 * @brief Strided index set start + stride * i, i in [0, size), as boost::numeric::ublas::slice
 * @details It stores the first index, the stride (a signed value) and the number of indices. The sentinel returned by
 * all() selects the whole extent of the target (with stride 1), and it is resolved by preprocess().
 */
class slice
{
public:
    /// Type of the indices and sizes of the slice. Alias of std::size_t.
    using size_type = std::size_t;

    /// Type of the stride of the slice. Alias of std::ptrdiff_t, the signed type of the stride.
    using difference_type = std::ptrdiff_t;

    /**
     * @brief Constructs the slice from its first index, its stride and its number of indices
     * @details The i-th element of the slice addresses the index Start + Stride * i.
     * @param Start First index of the slice.
     * @param Stride Signed step between two consecutive indices of the slice.
     * @param Size Number of indices in the slice.
     */
    slice(const size_type Start, const difference_type Stride, const size_type Size) : mStart(Start), mStride(Stride), mSize(Size) {}

    /**
     * @brief First index of the slice
     * @details Returns the start given at construction.
     * @return First index of the slice.
     */
    size_type start() const { return mStart; }

    /**
     * @brief Stride between consecutive elements of the slice
     * @details Returns the stride given at construction.
     * @return Stride between two consecutive indices of the slice.
     */
    difference_type stride() const { return mStride; }

    /**
     * @brief Number of indices in the slice
     * @details Returns the size given at construction.
     * @return Number of indices in the slice.
     */
    size_type size() const { return mSize; }

    /**
     * @brief Checks whether the slice contains no index
     * @details Returns true when the size is zero.
     * @return True if the slice contains no index, false otherwise.
     */
    bool empty() const { return mSize == 0; }

    /**
     * @brief Index of the target addressed by the Index-th element of the slice
     * @details Returns start + stride * Index, with the product computed on the signed stride and with no bound check.
     * @param Index Position inside the slice, counted from its first element.
     * @return Index of the target, that is start() + stride() * Index.
     */
    size_type operator()(const size_type Index) const { return mStart + static_cast<size_type>(mStride * static_cast<difference_type>(Index)); }

    /**
     * @brief Index of the target addressed by the Index-th element of the slice
     * @details Same as operator().
     * @param Index Position inside the slice, counted from its first element.
     * @return Index of the target, that is start() + stride() * Index.
     */
    size_type operator[](const size_type Index) const { return (*this)(Index); }

    /**
     * @brief Sentinel selecting the whole extent of the target
     * @details It is the slice starting at 0 with stride 1 and the maximum representable size, and it is resolved
     * against the actual extent by preprocess().
     * @return Sentinel slice starting at 0 with stride 1 and the maximum representable size.
     */
    static slice all() { return slice(0, 1, std::numeric_limits<size_type>::max()); }

    /**
     * @brief Slice resolved against the extent of the target it is applied to
     * @details Returns slice(0, 1, Extent) when this slice is the all() sentinel, otherwise this slice unchanged.
     * @param Extent Actual extent of the target the slice is applied to.
     * @return slice(0, 1, Extent) when this slice is the all() sentinel, otherwise this slice unchanged.
     */
    slice preprocess(const size_type Extent) const
    {
        return (mStart == 0 && mStride == 1 && mSize == std::numeric_limits<size_type>::max()) ? slice(0, 1, Extent) : *this;
    }

private:
    /**
     * @brief First index of the slice
     * @details Set at construction.
     */
    size_type mStart;

    /**
     * @brief Stride between consecutive elements of the slice
     * @details Set at construction.
     */
    difference_type mStride;

    /**
     * @brief Number of indices in the slice
     * @details Set at construction.
     */
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
/// Eigen expression type underlying matrix_row. The transposed row block of TExpression (a column-shaped expression),
/// with the cv-qualifiers removed.
template<class TExpression> using MatrixRowBase = std::remove_cv_t<decltype(std::declval<TExpression&>().row(0).transpose())>;

/// Eigen expression type underlying matrix_column. The column block of TExpression, with the cv-qualifiers removed.
template<class TExpression> using MatrixColumnBase = std::remove_cv_t<decltype(std::declval<TExpression&>().col(0))>;

/// Eigen expression type underlying vector_range. The segment of TExpression, with the cv-qualifiers removed.
template<class TExpression> using VectorRangeBase = std::remove_cv_t<decltype(std::declval<TExpression&>().segment(0, 0))>;

/// Eigen expression type underlying matrix_range. The block of TExpression, with the cv-qualifiers removed.
template<class TExpression> using MatrixRangeBase = std::remove_cv_t<decltype(std::declval<TExpression&>().block(0, 0, 0, 0))>;

/// Eigen expression type underlying vector_slice. An unaligned Eigen::Map with a dynamic inner stride over a dynamic
/// column vector of the value type of TExpression, the mapped vector being const when TExpression is const.
template<class TExpression> using VectorSliceBase = Eigen::Map<
    std::conditional_t<std::is_const_v<TExpression>, const DynamicColumnVector<typename TExpression::value_type>, DynamicColumnVector<typename TExpression::value_type>>,
    Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>;
} // namespace Internals

/**
 * @brief Row proxy (boost::numeric::ublas::matrix_row): a column-shaped, writable view of one matrix row
 * @details It derives from the Eigen expression of the transposed row, so it references the storage of the matrix
 * (nothing is copied) and every Eigen operation applies to it.
 * @tparam TExpression Type of the matrix, const-qualified for a read-only view.
 */
template<class TExpression>
class matrix_row : public Internals::MatrixRowBase<TExpression>
{
public:
    /// Eigen expression type this proxy derives from. The transposed row block of TExpression, see
    /// Internals::MatrixRowBase.
    using BaseType = Internals::MatrixRowBase<TExpression>;

    /// Scalar type of the entries of the view. Taken from the Eigen base expression.
    using value_type = typename BaseType::Scalar;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the view of the row Index of a matrix
     * @details The view is the transpose of the row of rMatrix, so it references the storage of the matrix and nothing
     * is copied.
     * @param rMatrix Matrix whose row is viewed.
     * @param Index Index of the row to view.
     */
    matrix_row(TExpression& rMatrix, const std::size_t Index) : BaseType(rMatrix.row(static_cast<Eigen::Index>(Index)).transpose()) {}

    /**
     * @brief Constructs the view from the underlying Eigen expression
     * @details Only the lightweight expression is copied, the entries stay in the storage it refers to.
     * @param rExpression Underlying Eigen expression, that is the transposed row block, to copy.
     */
    matrix_row(const BaseType& rExpression) : BaseType(rExpression) {}

    /**
     * @brief Copy constructor
     * @details The copy refers to the same row of the same matrix.
     * @param rOther Row proxy to copy.
     */
    matrix_row(const matrix_row& rOther) = default;

    /**
     * @brief Copy assignment from another row proxy
     * @details Copies the entries of rOther into the viewed row, the view itself is not rebound. The assignment from
     * any other Eigen expression is inherited from the Eigen base class by the using-declaration below.
     * @param rOther Row proxy whose entries are copied into the viewed row.
     * @return Reference to this view.
     */
    matrix_row& operator=(const matrix_row& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /**
     * @brief Assigns an expression to the view (uBLAS assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and assigned through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression assigned to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_row& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Adds an expression to the view (uBLAS plus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and added through Eigen's noalias(), so no temporary is created. Returns the
     * view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression added to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_row& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Subtracts an expression from the view (uBLAS minus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and subtracted through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression subtracted from the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_row& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Number of entries of the view
     * @details Returns the size of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of entries of the view.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/**
 * @brief Column proxy (boost::numeric::ublas::matrix_column): a writable view of one matrix column
 * @details It derives from the Eigen column block, so it references the storage of the matrix (nothing is copied) and
 * every Eigen operation applies to it.
 * @tparam TExpression Type of the matrix, const-qualified for a read-only view.
 */
template<class TExpression>
class matrix_column : public Internals::MatrixColumnBase<TExpression>
{
public:
    /// Eigen expression type this proxy derives from. The column block of TExpression, see Internals::MatrixColumnBase.
    using BaseType = Internals::MatrixColumnBase<TExpression>;

    /// Scalar type of the entries of the view. Taken from the Eigen base expression.
    using value_type = typename BaseType::Scalar;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the view of the column Index of a matrix
     * @details The view is the column of rMatrix, so it references the storage of the matrix and nothing is copied.
     * @param rMatrix Matrix whose column is viewed.
     * @param Index Index of the column to view.
     */
    matrix_column(TExpression& rMatrix, const std::size_t Index) : BaseType(rMatrix.col(static_cast<Eigen::Index>(Index))) {}

    /**
     * @brief Constructs the view from the underlying Eigen expression
     * @details Only the lightweight expression is copied, the entries stay in the storage it refers to.
     * @param rExpression Underlying Eigen expression, that is the column block, to copy.
     */
    matrix_column(const BaseType& rExpression) : BaseType(rExpression) {}

    /**
     * @brief Copy constructor
     * @details The copy refers to the same column of the same matrix.
     * @param rOther Column proxy to copy.
     */
    matrix_column(const matrix_column& rOther) = default;

    /**
     * @brief Copy assignment from another column proxy
     * @details Copies the entries of rOther into the viewed column, the view itself is not rebound. The assignment from
     * any other Eigen expression is inherited from the Eigen base class by the using-declaration below.
     * @param rOther Column proxy whose entries are copied into the viewed column.
     * @return Reference to this view.
     */
    matrix_column& operator=(const matrix_column& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /**
     * @brief Assigns an expression to the view (uBLAS assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and assigned through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression assigned to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_column& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Adds an expression to the view (uBLAS plus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and added through Eigen's noalias(), so no temporary is created. Returns the
     * view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression added to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_column& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Subtracts an expression from the view (uBLAS minus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and subtracted through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression subtracted from the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_column& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Number of entries of the view
     * @details Returns the size of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of entries of the view.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/**
 * @brief Range proxy (boost::numeric::ublas::vector_range): a writable view of a contiguous vector subrange
 * @details It derives from the Eigen segment of the vector, so it references the storage of the vector (nothing is
 * copied) and every Eigen operation applies to it.
 * @tparam TExpression Type of the vector, const-qualified for a read-only view.
 */
template<class TExpression>
class vector_range : public Internals::VectorRangeBase<TExpression>
{
public:
    /// Eigen expression type this proxy derives from. The segment of TExpression, see Internals::VectorRangeBase.
    using BaseType = Internals::VectorRangeBase<TExpression>;

    /// Scalar type of the entries of the view. Taken from the Eigen base expression.
    using value_type = typename BaseType::Scalar;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the view of the entries of a vector selected by a range
     * @details The range is resolved against the size of the vector (range::all() selects the whole vector). The view
     * references the storage of the vector and nothing is copied.
     * @param rVector Vector whose entries are viewed.
     * @param rRange Range of entries to view, resolved against the size of the vector.
     */
    vector_range(TExpression& rVector, const range& rRange)
        : BaseType(rVector.segment(static_cast<Eigen::Index>(rRange.preprocess(rVector.size()).start()), static_cast<Eigen::Index>(rRange.preprocess(rVector.size()).size()))) {}

    /**
     * @brief Constructs the view of the entries [Start, Stop) of a vector
     * @details Delegates to the constructor taking the range(Start, Stop).
     * @param rVector Vector whose entries are viewed.
     * @param Start First index of the view.
     * @param Stop End of the view, the first index that does not belong to it.
     */
    vector_range(TExpression& rVector, const std::size_t Start, const std::size_t Stop) : vector_range(rVector, range(Start, Stop)) {}

    /**
     * @brief Constructs the view from the underlying Eigen expression
     * @details Only the lightweight expression is copied, the entries stay in the storage it refers to.
     * @param rExpression Underlying Eigen expression, that is the segment, to copy.
     */
    vector_range(const BaseType& rExpression) : BaseType(rExpression) {}

    /**
     * @brief Copy constructor
     * @details The copy refers to the same subrange of the same vector.
     * @param rOther Range proxy to copy.
     */
    vector_range(const vector_range& rOther) = default;

    /**
     * @brief Copy assignment from another range proxy
     * @details Copies the entries of rOther into the viewed subrange, the view itself is not rebound. The assignment
     * from any other Eigen expression is inherited from the Eigen base class by the using-declaration below.
     * @param rOther Range proxy whose entries are copied into the viewed subrange.
     * @return Reference to this view.
     */
    vector_range& operator=(const vector_range& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /**
     * @brief Assigns an expression to the view (uBLAS assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and assigned through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression assigned to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_range& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Adds an expression to the view (uBLAS plus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and added through Eigen's noalias(), so no temporary is created. Returns the
     * view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression added to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_range& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Subtracts an expression from the view (uBLAS minus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and subtracted through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression subtracted from the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_range& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Number of entries of the view
     * @details Returns the size of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of entries of the view.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/**
 * @brief Range proxy (boost::numeric::ublas::matrix_range): a writable view of a matrix block
 * @details It derives from the Eigen block, so it references the storage of the matrix (nothing is copied) and every
 * Eigen operation applies to it.
 * @tparam TExpression Type of the matrix, const-qualified for a read-only view.
 */
template<class TExpression>
class matrix_range : public Internals::MatrixRangeBase<TExpression>
{
public:
    /// Eigen expression type this proxy derives from. The block of TExpression, see Internals::MatrixRangeBase.
    using BaseType = Internals::MatrixRangeBase<TExpression>;

    /// Scalar type of the entries of the view. Taken from the Eigen base expression.
    using value_type = typename BaseType::Scalar;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the view of the block of a matrix selected by two ranges
     * @details The first range selects the rows and the second one the columns, each resolved against the corresponding
     * extent of the matrix (range::all() selects the whole extent). The view references the storage of the matrix and
     * nothing is copied.
     * @param rMatrix Matrix whose block is viewed.
     * @param rRange1 Range of rows to view, resolved against the number of rows of the matrix.
     * @param rRange2 Range of columns to view, resolved against the number of columns of the matrix.
     */
    matrix_range(TExpression& rMatrix, const range& rRange1, const range& rRange2)
        : BaseType(rMatrix.block(
            static_cast<Eigen::Index>(rRange1.preprocess(rMatrix.rows()).start()),
            static_cast<Eigen::Index>(rRange2.preprocess(rMatrix.cols()).start()),
            static_cast<Eigen::Index>(rRange1.preprocess(rMatrix.rows()).size()),
            static_cast<Eigen::Index>(rRange2.preprocess(rMatrix.cols()).size()))) {}

    /**
     * @brief Constructs the view from the underlying Eigen expression
     * @details Only the lightweight expression is copied, the entries stay in the storage it refers to.
     * @param rExpression Underlying Eigen expression, that is the block, to copy.
     */
    matrix_range(const BaseType& rExpression) : BaseType(rExpression) {}

    /**
     * @brief Copy constructor
     * @details The copy refers to the same block of the same matrix.
     * @param rOther Range proxy to copy.
     */
    matrix_range(const matrix_range& rOther) = default;

    /**
     * @brief Copy assignment from another range proxy
     * @details Copies the entries of rOther into the viewed block, the view itself is not rebound. The assignment from
     * any other Eigen expression is inherited from the Eigen base class by the using-declaration below.
     * @param rOther Range proxy whose entries are copied into the viewed block.
     * @return Reference to this view.
     */
    matrix_range& operator=(const matrix_range& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /**
     * @brief Assigns an expression to the view (uBLAS assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and assigned through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression assigned to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_range& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Adds an expression to the view (uBLAS plus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and added through Eigen's noalias(), so no temporary is created. Returns the
     * view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression added to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_range& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Subtracts an expression from the view (uBLAS minus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and subtracted through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression subtracted from the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> matrix_range& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Number of rows of the view
     * @details Returns the number of rows of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of rows of the view.
     */
    std::size_t size1() const { return static_cast<std::size_t>(BaseType::rows()); }

    /**
     * @brief Number of columns of the view
     * @details Returns the number of columns of the Eigen base expression as an unsigned integer, as in the uBLAS
     * surface.
     * @return Number of columns of the view.
     */
    std::size_t size2() const { return static_cast<std::size_t>(BaseType::cols()); }
};

/**
 * @brief Slice proxy (boost::numeric::ublas::vector_slice): a writable strided view of a vector
 * @details It derives from an Eigen::Map with a dynamic inner stride over the contiguous storage of the vector, so it
 * references that storage (nothing is copied) and every Eigen operation applies to it.
 * @tparam TExpression Type of the vector, const-qualified for a read-only view.
 */
template<class TExpression>
class vector_slice : public Internals::VectorSliceBase<TExpression>
{
public:
    /// Eigen expression type this proxy derives from. An Eigen::Map with a dynamic inner stride over the storage of the
    /// vector, see Internals::VectorSliceBase.
    using BaseType = Internals::VectorSliceBase<TExpression>;

    /// Scalar type of the entries of the view. Taken from the Eigen base expression.
    using value_type = typename BaseType::Scalar;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the strided view of the entries of a vector selected by a slice
     * @details The slice is resolved against the size of the vector (slice::all() selects the whole vector). The view
     * is a map over the storage of the vector that starts at the first index of the slice and uses its stride as inner
     * stride, so nothing is copied.
     * @param rVector Vector whose entries are viewed.
     * @param rSlice Slice of entries to view, resolved against the size of the vector.
     */
    vector_slice(TExpression& rVector, const slice& rSlice)
        : BaseType(rVector.data().begin() + rSlice.preprocess(rVector.size()).start(),
                   static_cast<Eigen::Index>(rSlice.preprocess(rVector.size()).size()),
                   Eigen::InnerStride<Eigen::Dynamic>(static_cast<Eigen::Index>(rSlice.preprocess(rVector.size()).stride()))) {}

    /**
     * @brief Constructs the view from the underlying Eigen map
     * @details Only the lightweight map is copied, the entries stay in the storage it refers to.
     * @param rExpression Underlying Eigen map to copy.
     */
    vector_slice(const BaseType& rExpression) : BaseType(rExpression) {}

    /**
     * @brief Copy constructor
     * @details The copy refers to the same slice of the same vector.
     * @param rOther Slice proxy to copy.
     */
    vector_slice(const vector_slice& rOther) = default;

    /**
     * @brief Copy assignment from another slice proxy
     * @details Copies the entries of rOther into the viewed slice, the view itself is not rebound. The assignment from
     * any other Eigen expression is inherited from the Eigen base class by the using-declaration below.
     * @param rOther Slice proxy whose entries are copied into the viewed slice.
     * @return Reference to this view.
     */
    vector_slice& operator=(const vector_slice& rOther) { BaseType::operator=(rOther); return *this; }
    using BaseType::operator=;

    /**
     * @brief Assigns an expression to the view (uBLAS assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and assigned through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression assigned to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_slice& assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() = rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Adds an expression to the view (uBLAS plus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and added through Eigen's noalias(), so no temporary is created. Returns the
     * view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression added to the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_slice& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() += rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Subtracts an expression from the view (uBLAS minus_assign)
     * @details Part of the uBLAS assignment protocol (assign/plus_assign/minus_assign) on the view. The expression is
     * cast to the scalar type of the view and subtracted through Eigen's noalias(), so no temporary is created. Returns
     * the view.
     * @tparam TDerived Derived type of the Eigen expression.
     * @param rExpression Expression subtracted from the view, cast to the scalar type of the view.
     * @return Reference to this view.
     */
    template<class TDerived> vector_slice& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression) { this->noalias() -= rExpression.derived().template cast<value_type>(); return *this; }

    /**
     * @brief Number of entries of the view
     * @details Returns the size of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of entries of the view.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }
};

/**
 * @brief Row of a matrix as a (column-shaped) vector proxy, readable and writable (ublas row())
 * @details Returns a matrix_row referencing the row I of rM, so nothing is copied. Only enabled for dense Eigen types
 * (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose row is viewed.
 * @param I Index of the row.
 * @return Writable proxy referencing the row I of rM.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_row<TMatrix> row(TMatrix& rM, const std::size_t I) { return matrix_row<TMatrix>(rM, I); }

/**
 * @brief Row of a matrix as a (column-shaped) vector proxy, read-only (ublas row())
 * @details Returns a matrix_row over the const matrix referencing the row I of rM, so nothing is copied. Only enabled
 * for dense Eigen types (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose row is viewed.
 * @param I Index of the row.
 * @return Read-only proxy referencing the row I of rM.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_row<const TMatrix> row(const TMatrix& rM, const std::size_t I) { return matrix_row<const TMatrix>(rM, I); }

/**
 * @brief Row of a compressed matrix (ublas matrix_row over a compressed_matrix)
 * @details Covers the uBLAS idioms used on sparse rows: zeroing the stored entries of the row (row(A, i) =
 * ZeroVector(n)) and the scalar product with a dense vector (inner_prod(v, row(A, i))). The sparsity pattern is kept.
 * @tparam TMatrix Compressed matrix type, possibly const-qualified.
 */
template<class TMatrix>
class compressed_matrix_row
{
public:
    /// Scalar type of the entries of the matrix. Taken from the matrix type without its const qualifier.
    using value_type = typename std::remove_const_t<TMatrix>::value_type;

    /**
     * @brief Constructs the proxy of the row Index of a compressed matrix
     * @details The matrix is referenced, not copied, so it must outlive the proxy.
     * @param rMatrix Compressed matrix whose row is viewed.
     * @param Index Index of the row.
     */
    compressed_matrix_row(TMatrix& rMatrix, const std::size_t Index) : mrMatrix(rMatrix), mIndex(Index) {}

    /**
     * @brief Length of the row
     * @details Returns the number of columns of the matrix.
     * @return Number of columns of the matrix, that is the length of the row.
     */
    std::size_t size() const { return mrMatrix.size2(); }

    /**
     * @brief Zeroes the stored entries of the row
     * @details Implements the uBLAS idiom row(A, i) = ZeroVector(n): the argument only selects this overload and its
     * size is not used. The stored values of the row are set to zero and the sparsity pattern is kept.
     * @return Reference to this proxy.
     */
    compressed_matrix_row& operator=(const zero_vector<value_type>&)
    {
        value_type* p_values = mrMatrix.valuePtr();
        for (std::size_t k = RowBegin(); k < RowEnd(); ++k) {
            p_values[k] = value_type();
        }
        return *this;
    }

    /**
     * @brief Scalar product with a dense vector expression, over the stored entries
     * @details Sums the products of each stored entry of the row by the entry of the vector at the column of that
     * entry. It backs the inner_prod() overloads taking a compressed_matrix_row.
     * @tparam TDerived Derived type of the dense Eigen vector expression.
     * @param rVector Dense vector expression, read at the column of each stored entry of the row.
     * @return Scalar product of the row with rVector, computed over the stored entries of the row.
     */
    template<class TDerived>
    value_type Dot(const Eigen::MatrixBase<TDerived>& rVector) const
    {
        const value_type* p_values = mrMatrix.valuePtr();
        const typename std::remove_const_t<TMatrix>::StorageIndex* p_columns = mrMatrix.innerIndexPtr();
        value_type result = value_type();
        for (std::size_t k = RowBegin(); k < RowEnd(); ++k) {
            result += p_values[k] * rVector.derived().coeff(static_cast<Eigen::Index>(p_columns[k]));
        }
        return result;
    }

private:
    /**
     * @brief Position of the first stored entry of the row
     * @details Position in the compressed arrays of values and column indices, read from the outer index array of the
     * matrix.
     * @return Position of the first stored entry of the row in the compressed arrays.
     */
    std::size_t RowBegin() const { return static_cast<std::size_t>(mrMatrix.outerIndexPtr()[mIndex]); }

    /**
     * @brief Position past the last stored entry of the row
     * @details Position in the compressed arrays of values and column indices, read from the outer index array of the
     * matrix.
     * @return Position past the last stored entry of the row in the compressed arrays.
     */
    std::size_t RowEnd() const { return static_cast<std::size_t>(mrMatrix.outerIndexPtr()[mIndex + 1]); }

    /**
     * @brief Compressed matrix the row belongs to
     * @details Referenced, not copied.
     */
    TMatrix& mrMatrix;

    /**
     * @brief Index of the row in the matrix
     * @details Set at construction.
     */
    std::size_t mIndex;
};

/**
 * @brief Row of a compressed matrix, zeroable and usable in inner_prod (ublas row())
 * @details Returns a compressed_matrix_row referencing the row I of rM. The proxy can be zeroed (row(A, i) =
 * ZeroVector(n)) and used in inner_prod().
 * @tparam TDataType Scalar type of the entries of the compressed matrix.
 * @tparam TIndexType Index type of the storage of the compressed matrix.
 * @param rM Compressed matrix whose row is viewed.
 * @param I Index of the row.
 * @return Proxy referencing the row I of rM, which can be zeroed and used in inner_prod().
 */
template<class TDataType, class TIndexType>
inline compressed_matrix_row<EigenCompressedMatrix<TDataType, TIndexType>> row(EigenCompressedMatrix<TDataType, TIndexType>& rM, const std::size_t I)
{
    return compressed_matrix_row<EigenCompressedMatrix<TDataType, TIndexType>>(rM, I);
}

/**
 * @brief Row of a compressed matrix, read-only (ublas row())
 * @details Returns a compressed_matrix_row over the const matrix referencing the row I of rM. The proxy can be used in
 * inner_prod() but its entries cannot be modified.
 * @tparam TDataType Scalar type of the entries of the compressed matrix.
 * @tparam TIndexType Index type of the storage of the compressed matrix.
 * @param rM Compressed matrix whose row is viewed.
 * @param I Index of the row.
 * @return Read-only proxy referencing the row I of rM, usable in inner_prod().
 */
template<class TDataType, class TIndexType>
inline compressed_matrix_row<const EigenCompressedMatrix<TDataType, TIndexType>> row(const EigenCompressedMatrix<TDataType, TIndexType>& rM, const std::size_t I)
{
    return compressed_matrix_row<const EigenCompressedMatrix<TDataType, TIndexType>>(rM, I);
}

/**
 * @brief Column of a matrix as a vector proxy, readable and writable (ublas column())
 * @details Returns a matrix_column referencing the column J of rM, so nothing is copied. Only enabled for dense Eigen
 * types (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose column is viewed.
 * @param J Index of the column.
 * @return Writable proxy referencing the column J of rM.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_column<TMatrix> column(TMatrix& rM, const std::size_t J) { return matrix_column<TMatrix>(rM, J); }

/**
 * @brief Column of a matrix as a vector proxy, read-only (ublas column())
 * @details Returns a matrix_column over the const matrix referencing the column J of rM, so nothing is copied. Only
 * enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose column is viewed.
 * @param J Index of the column.
 * @return Read-only proxy referencing the column J of rM.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_column<const TMatrix> column(const TMatrix& rM, const std::size_t J) { return matrix_column<const TMatrix>(rM, J); }

/**
 * @brief Vector subrange [Low, High), readable and writable (ublas subrange())
 * @details Returns a vector_range referencing the entries of rV from Low up to High (excluded), so nothing is copied.
 * Only enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TVector Dense Eigen vector type.
 * @param rV Vector whose entries are viewed.
 * @param Low First index of the subrange.
 * @param High End of the subrange, the first index that does not belong to it.
 * @return Writable proxy referencing the entries of rV from Low up to High (excluded).
 */
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<TVector> subrange(TVector& rV, const std::size_t Low, const std::size_t High) { return vector_range<TVector>(rV, Low, High); }

/**
 * @brief Vector subrange [Low, High), read-only (ublas subrange())
 * @details Returns a vector_range over the const vector referencing the entries of rV from Low up to High (excluded),
 * so nothing is copied. Only enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TVector Dense Eigen vector type.
 * @param rV Vector whose entries are viewed.
 * @param Low First index of the subrange.
 * @param High End of the subrange, the first index that does not belong to it.
 * @return Read-only proxy referencing the entries of rV from Low up to High (excluded).
 */
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<const TVector> subrange(const TVector& rV, const std::size_t Low, const std::size_t High) { return vector_range<const TVector>(rV, Low, High); }

/**
 * @brief Matrix block [Row1, Row2) x [Col1, Col2), readable and writable (ublas subrange())
 * @details Returns a matrix_range referencing the rows of rM from Row1 up to Row2 (excluded) and its columns from Col1
 * up to Col2 (excluded), so nothing is copied. Only enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose block is viewed.
 * @param Row1 First row of the block.
 * @param Row2 End of the rows of the block, the first row that does not belong to it.
 * @param Col1 First column of the block.
 * @param Col2 End of the columns of the block, the first column that does not belong to it.
 * @return Writable proxy referencing the block of rM between the rows Row1 and Row2 (excluded) and the columns Col1 and
 * Col2 (excluded).
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<TMatrix> subrange(TMatrix& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return matrix_range<TMatrix>(rM, range(Row1, Row2), range(Col1, Col2));
}

/**
 * @brief Matrix block [Row1, Row2) x [Col1, Col2), read-only (ublas subrange())
 * @details Returns a matrix_range over the const matrix referencing the rows of rM from Row1 up to Row2 (excluded) and
 * its columns from Col1 up to Col2 (excluded), so nothing is copied. Only enabled for dense Eigen types
 * (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose block is viewed.
 * @param Row1 First row of the block.
 * @param Row2 End of the rows of the block, the first row that does not belong to it.
 * @param Col1 First column of the block.
 * @param Col2 End of the columns of the block, the first column that does not belong to it.
 * @return Read-only proxy referencing the block of rM between the rows Row1 and Row2 (excluded) and the columns Col1
 * and Col2 (excluded).
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<const TMatrix> subrange(const TMatrix& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return matrix_range<const TMatrix>(rM, range(Row1, Row2), range(Col1, Col2));
}

/**
 * @brief Vector range projection, readable and writable (ublas project())
 * @details Returns a vector_range referencing the entries of rV selected by rRange, which is resolved against the size
 * of the vector, so nothing is copied. Only enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TVector Dense Eigen vector type.
 * @param rV Vector whose entries are viewed.
 * @param rRange Range of entries to select, resolved against the size of the vector.
 * @return Writable proxy referencing the entries of rV selected by rRange.
 */
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<TVector> project(TVector& rV, const range& rRange) { return vector_range<TVector>(rV, rRange); }

/**
 * @brief Vector range projection, read-only (ublas project())
 * @details Returns a vector_range over the const vector referencing the entries of rV selected by rRange, which is
 * resolved against the size of the vector, so nothing is copied. Only enabled for dense Eigen types
 * (Internals::IsEigenDense).
 * @tparam TVector Dense Eigen vector type.
 * @param rV Vector whose entries are viewed.
 * @param rRange Range of entries to select, resolved against the size of the vector.
 * @return Read-only proxy referencing the entries of rV selected by rRange.
 */
template<class TVector> requires Internals::IsEigenDense<TVector>
inline vector_range<const TVector> project(const TVector& rV, const range& rRange) { return vector_range<const TVector>(rV, rRange); }

/**
 * @brief Matrix range projection, readable and writable (ublas project())
 * @details Returns a matrix_range referencing the block of rM whose rows are selected by rRange1 and whose columns are
 * selected by rRange2, so nothing is copied. Only enabled for dense Eigen types (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose block is viewed.
 * @param rRange1 Range of rows to select, resolved against the number of rows of the matrix.
 * @param rRange2 Range of columns to select, resolved against the number of columns of the matrix.
 * @return Writable proxy referencing the block of rM selected by rRange1 and rRange2.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<TMatrix> project(TMatrix& rM, const range& rRange1, const range& rRange2) { return matrix_range<TMatrix>(rM, rRange1, rRange2); }

/**
 * @brief Matrix range projection, read-only (ublas project())
 * @details Returns a matrix_range over the const matrix referencing the block of rM whose rows are selected by rRange1
 * and whose columns are selected by rRange2, so nothing is copied. Only enabled for dense Eigen types
 * (Internals::IsEigenDense).
 * @tparam TMatrix Dense Eigen matrix type.
 * @param rM Matrix whose block is viewed.
 * @param rRange1 Range of rows to select, resolved against the number of rows of the matrix.
 * @param rRange2 Range of columns to select, resolved against the number of columns of the matrix.
 * @return Read-only proxy referencing the block of rM selected by rRange1 and rRange2.
 */
template<class TMatrix> requires Internals::IsEigenDense<TMatrix>
inline matrix_range<const TMatrix> project(const TMatrix& rM, const range& rRange1, const range& rRange2) { return matrix_range<const TMatrix>(rM, rRange1, rRange2); }

/**
 * @brief Vector slice projection over a storage-owning vector, readable and writable (ublas project())
 * @details Returns a vector_slice referencing the entries of rV selected by rSlice, so nothing is copied. Only enabled
 * for dense Eigen vectors whose base is a plain (storage-owning) Eigen object, since the slice is a strided map over
 * the contiguous storage.
 * @tparam TVector Dense Eigen vector type whose base is a plain (storage-owning) Eigen object.
 * @param rV Vector whose entries are viewed.
 * @param rSlice Slice of entries to select, resolved against the size of the vector.
 * @return Writable strided proxy referencing the entries of rV selected by rSlice.
 */
template<class TVector> requires (Internals::IsEigenDense<TVector> && Internals::IsPlainObject<typename TVector::BaseType>)
inline vector_slice<TVector> project(TVector& rV, const slice& rSlice) { return vector_slice<TVector>(rV, rSlice); }

/**
 * @brief Vector slice projection over a storage-owning vector, read-only (ublas project())
 * @details Returns a vector_slice over the const vector referencing the entries of rV selected by rSlice, so nothing is
 * copied. Only enabled for dense Eigen vectors whose base is a plain (storage-owning) Eigen object, since the slice is
 * a strided map over the contiguous storage.
 * @tparam TVector Dense Eigen vector type whose base is a plain (storage-owning) Eigen object.
 * @param rV Vector whose entries are viewed.
 * @param rSlice Slice of entries to select, resolved against the size of the vector.
 * @return Read-only strided proxy referencing the entries of rV selected by rSlice.
 */
template<class TVector> requires (Internals::IsEigenDense<TVector> && Internals::IsPlainObject<typename TVector::BaseType>)
inline vector_slice<const TVector> project(const TVector& rV, const slice& rSlice) { return vector_slice<const TVector>(rV, rSlice); }

///@}
///@name Products
///@{

namespace Internals
{
/**
 * @brief Shared implementation for the dense prod() overloads below
 * @details Kept out of the "prod" overload set on purpose: a call from inside
 * a prod() overload to unqualified prod(rA, rB) would re-run overload
 * resolution against every prod overload (including the sparse ones), and
 * deducing those against a reference-to-MatrixBase argument forces an
 * instantiation of Eigen::SparseMatrixBase<TDerived> that hard-errors.
 * A vector first operand keeps the ublas prod(v, M) semantics (v^T M),
 * returned column-shaped so it assigns to the vector types.
 * @tparam TDerived1 Derived type of the first operand.
 * @tparam TDerived2 Derived type of the second operand.
 * @param rA First operand of the product.
 * @param rB Second operand of the product.
 * @return Lazy Eigen expression of the product: the transpose of rA^T rB when rA is a vector, rA rB otherwise.
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
 * @brief Matrix/matrix and matrix/vector product (lazy Eigen expression)
 * @details The requires-clause guards TDerived1/TDerived2 the same way as the
 * sparse overloads below: prod<TResult>(...) shares the "prod" name, so a call
 * like prod<Matrix>(...) explicitly substitutes TDerived1 = Matrix into every
 * prod overload, and naming Eigen::MatrixBase<Matrix> for that substitution
 * forces Eigen to complete it, which hard-errors (Eigen::internal::traits is
 * unspecialized for the Kratos wrapper types).
 * @tparam TShield Leading parameter, never specified by callers, that absorbs the explicit first template argument of
 * the calls pinning the result type.
 * @tparam TDerived1 Derived type of the first operand.
 * @tparam TDerived2 Derived type of the second operand.
 * @param rA First operand of the product.
 * @param rB Second operand of the product.
 * @return Lazy Eigen expression of the product, see Internals::dense_prod().
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

/**
 * @brief Product with an explicitly pinned result type, as ublas prod<Vector>(A, B)
 * @details Evaluates the dense product (see the prod() overload above) and converts it to TResult, so the result is a
 * concrete object instead of a lazy Eigen expression.
 * @tparam TResult Concrete type the product is converted to.
 * @tparam TDerived1 Derived type of the first operand.
 * @tparam TDerived2 Derived type of the second operand.
 * @param rA First operand of the product.
 * @param rB Second operand of the product.
 * @return Product evaluated into an object of type TResult.
 */
template<class TResult, class TDerived1, class TDerived2>
inline TResult prod(const Eigen::MatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rB)
{
    return TResult(Internals::dense_prod(rA, rB));
}

/**
 * @brief Sparse-matrix/vector (or sparse/dense) product
 * @details Returns the lazy Eigen expression rA * rX. It is guarded like the dense prod() overload, see the TShield
 * parameter there.
 * @tparam TShield Leading guard parameter, never specified by callers, see the dense prod() overload.
 * @tparam TDerived1 Derived type of the sparse matrix.
 * @tparam TDerived2 Derived type of the dense vector or matrix.
 * @param rA Sparse matrix.
 * @param rX Dense vector or matrix multiplied by the sparse matrix.
 * @return Lazy Eigen expression rA * rX.
 */
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires { typename Eigen::internal::traits<TDerived1>::StorageKind; })
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::MatrixBase<TDerived2>& rX)
{
    return rA.derived() * rX.derived();
}

/**
 * @brief Vector/sparse-matrix product: v^T A, returned column-shaped (uBLAS prod(v, A))
 * @details Returns the lazy Eigen expression A^T v, the transpose of v^T A, so that the result is column-shaped and
 * assigns to the vector types. It is guarded like the dense prod() overload, see the TShield parameter there.
 * @tparam TShield Leading guard parameter, never specified by callers, see the dense prod() overload.
 * @tparam TDerived1 Derived type of the dense vector.
 * @tparam TDerived2 Derived type of the sparse matrix.
 * @param rX Dense vector multiplying the sparse matrix from the left.
 * @param rA Sparse matrix.
 * @return Lazy Eigen expression A^T v, the column-shaped transpose of v^T A.
 */
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires { typename Eigen::internal::traits<TDerived2>::StorageKind; })
inline auto prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::SparseMatrixBase<TDerived2>& rA)
{
    return rA.derived().transpose() * rX.derived();
}

/**
 * @brief Sparse-matrix/sparse-matrix product
 * @details Returns the lazy Eigen expression rA * rB. It is guarded like the dense prod() overload, see the TShield
 * parameter there.
 * @tparam TShield Leading guard parameter, never specified by callers, see the dense prod() overload.
 * @tparam TDerived1 Derived type of the first sparse matrix.
 * @tparam TDerived2 Derived type of the second sparse matrix.
 * @param rA First sparse matrix.
 * @param rB Second sparse matrix.
 * @return Lazy Eigen expression rA * rB.
 */
template<class TShield = void, class TDerived1, class TDerived2>
    requires (std::is_void_v<TShield> && requires {
        typename Eigen::internal::traits<TDerived1>::StorageKind;
        typename Eigen::internal::traits<TDerived2>::StorageKind;
    })
inline auto prod(const Eigen::SparseMatrixBase<TDerived1>& rA, const Eigen::SparseMatrixBase<TDerived2>& rB)
{
    return rA.derived() * rB.derived();
}

/**
 * @brief Scalar (dot) product of two vectors
 * @details As in uBLAS the product is NOT conjugated for complex scalars (Eigen's dot() would conjugate). It is
 * computed as the 1 x 1 product rX^T rY and, when the scalar types differ, both operands are first cast to the type of
 * their product.
 * @tparam TDerived1 Derived type of the first vector.
 * @tparam TDerived2 Derived type of the second vector.
 * @param rX First vector.
 * @param rY Second vector.
 * @return Scalar product of rX and rY, not conjugated for complex scalars.
 */
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

/**
 * @brief Scalar product of a dense vector and a row of a compressed matrix (ublas inner_prod)
 * @details Delegates to compressed_matrix_row::Dot(), so only the stored entries of the row contribute.
 * @tparam TDerived Derived type of the dense vector.
 * @tparam TMatrix Compressed matrix type of the row, possibly const-qualified.
 * @param rX Dense vector.
 * @param rY Row of a compressed matrix.
 * @return Scalar product of rX and the row rY, over the stored entries of the row.
 */
template<class TDerived, class TMatrix>
inline typename compressed_matrix_row<TMatrix>::value_type inner_prod(const Eigen::MatrixBase<TDerived>& rX, const compressed_matrix_row<TMatrix>& rY)
{
    return rY.Dot(rX);
}

/**
 * @brief Scalar product of a row of a compressed matrix and a dense vector (ublas inner_prod)
 * @details Delegates to compressed_matrix_row::Dot(), so only the stored entries of the row contribute.
 * @tparam TMatrix Compressed matrix type of the row, possibly const-qualified.
 * @tparam TDerived Derived type of the dense vector.
 * @param rX Row of a compressed matrix.
 * @param rY Dense vector.
 * @return Scalar product of the row rX and rY, over the stored entries of the row.
 */
template<class TMatrix, class TDerived>
inline typename compressed_matrix_row<TMatrix>::value_type inner_prod(const compressed_matrix_row<TMatrix>& rX, const Eigen::MatrixBase<TDerived>& rY)
{
    return rX.Dot(rY);
}

/**
 * @brief Outer product of two vectors, rX * rY^T
 * @details A row-shaped rY (as produced by Eigen's own transpose) is accepted as is, otherwise rY is transposed. The
 * result is a lazy Eigen expression.
 * @tparam TDerived1 Derived type of the first vector.
 * @tparam TDerived2 Derived type of the second vector.
 * @param rX First vector, the column factor of the product.
 * @param rY Second vector, transposed in the product unless it is already row-shaped.
 * @return Lazy Eigen expression of the matrix rX * rY^T.
 */
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
 * @brief uBLAS-style axpy product: rY = A B (Init = true) or rY += A B
 * @details A single overload covers the three uBLAS forms, dispatching on the
 * compile-time shape of the first operand exactly as prod() does: matrix x
 * vector, matrix x matrix and, for a vector first operand, the vector x
 * matrix form (y = M^T v). As in uBLAS the target must not alias either
 * operand; the assignment goes through Eigen's noalias() path.
 * @tparam TDerived1 Derived type of the first operand.
 * @tparam TDerived2 Derived type of the second operand.
 * @tparam TTargetType Type of the target, a dense Eigen object.
 * @param rA First operand of the product.
 * @param rB Second operand of the product.
 * @param rY Target of the product, it must not alias either operand.
 * @param Init If true the target is overwritten with the product, otherwise the product is added to it.
 * @return Reference to the target rY.
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

/**
 * @brief Sparse counterpart of the axpy product (sparse matrix times a dense vector or matrix)
 * @details Computes rY = A X (Init = true) or rY += A X through Eigen's noalias() path and returns rY. As for the dense
 * overload, the target must not alias the operands.
 * @tparam TDerived1 Derived type of the sparse matrix.
 * @tparam TDerived2 Derived type of the dense vector or matrix.
 * @tparam TTargetType Type of the target, a dense Eigen object.
 * @param rA Sparse matrix.
 * @param rX Dense vector or matrix multiplied by the sparse matrix.
 * @param rY Target of the product, it must not alias either operand.
 * @param Init If true the target is overwritten with the product, otherwise the product is added to it.
 * @return Reference to the target rY.
 */
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
 * @brief Transpose, with the uBLAS TYPE-based semantics
 * @details uBLAS defines trans() on a vector as the identity ((trans v)[i] = v[i],
 * boost vector_expression.hpp) and discriminates vector from matrix by type,
 * not by shape. Every compile-time vector-shaped expression (Vector, array_1d,
 * BoundedVector, the row/column proxies, ...) therefore comes back unchanged
 * (as a const reference to the argument, so no copy is made), so idioms like outer_prod(x, trans(y)) and
 * prod(trans(v), M) mean the same under both backends; matrices are lazily transposed.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Expression to transpose.
 * @return Const reference to the argument when it is column-shaped, the lazy transpose of the argument otherwise.
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

/**
 * @brief Transpose of a bounded matrix, whatever its static shape
 * @details A bounded matrix is a matrix by type: an N x 1 or 1 x N BoundedMatrix transposes (unlike the N-vectors, see
 * the generic overload above). This overload is an exact match for the wrapper type, so it beats the generic MatrixBase
 * template.
 * @tparam T Scalar type of the bounded matrix.
 * @tparam TSize1 Fixed number of rows of the bounded matrix.
 * @tparam TSize2 Fixed number of columns of the bounded matrix.
 * @param rM Bounded matrix to transpose.
 * @return Lazy transpose of the matrix.
 */
template<class T, std::size_t TSize1, std::size_t TSize2>
inline auto trans(const EigenBoundedMatrix<T, TSize1, TSize2>& rM)
{
    return rM.transpose();
}

/**
 * @brief Transpose of a sparse matrix
 * @details Unlike the dense overload this MATERIALIZES the result: Eigen's lazy sparse transpose flips the storage
 * order, and sparse binary operations require both sides to share the same order.
 * @tparam TDerived Derived type of the sparse matrix.
 * @param rM Sparse matrix to transpose.
 * @return Materialized transpose of the matrix, as a plain object.
 */
template<class TDerived>
inline typename TDerived::PlainObject trans(const Eigen::SparseMatrixBase<TDerived>& rM)
{
    return rM.transpose();
}

///@}
///@name Norms and reductions
///@{

/**
 * @brief 1-norm: the sum of absolute values for a vector, the maximum absolute column sum for a matrix
 * @details ublas semantics: an expression that is a column at compile time (Internals::IsColumnShaped) is a vector, any
 * other one is a matrix. The result is 0 for an empty matrix.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Vector or matrix.
 * @return 1-norm of rX: the sum of the absolute values for a vector, the maximum absolute column sum for a matrix.
 */
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

/**
 * @brief Euclidean norm (Frobenius norm for matrices, as in ublas)
 * @details Returns the norm() of the Eigen expression, which is the Frobenius norm for a matrix.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Vector or matrix.
 * @return Euclidean norm of rX, which is the Frobenius norm for a matrix.
 */
template<class TDerived>
inline auto norm_2(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.norm();
}

/**
 * @brief Squared Euclidean norm
 * @details Returns the squaredNorm() of the Eigen expression, which is the squared Frobenius norm for a matrix.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Vector or matrix.
 * @return Squared Euclidean norm of rX, which is the squared Frobenius norm for a matrix.
 */
template<class TDerived>
inline auto norm_2_square(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.squaredNorm();
}

/**
 * @brief Infinity norm: the maximum absolute value for a vector, the maximum absolute row sum for a matrix
 * @details ublas semantics: an expression that is a column at compile time (Internals::IsColumnShaped) is a vector, any
 * other one is a matrix. The result is 0 for an empty matrix.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Vector or matrix.
 * @return Infinity norm of rX: the maximum absolute value for a vector, the maximum absolute row sum for a matrix.
 */
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

/**
 * @brief Index of the first entry of maximum absolute value (ublas index_norm_inf)
 * @details The entries are visited with linear indexing and, in case of ties, the first one is kept. The first entry is
 * read unconditionally, so the argument is expected to be a non-empty vector.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Non-empty vector.
 * @return Index of the first entry of maximum absolute value.
 */
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

/**
 * @brief Sum of all coefficients
 * @details Returns the sum() of the Eigen expression.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rX Vector or matrix.
 * @return Sum of all the coefficients of rX.
 */
template<class TDerived>
inline auto sum(const Eigen::MatrixBase<TDerived>& rX)
{
    return rX.sum();
}

/**
 * @brief Frobenius norm (same as norm_2 for vectors, as in ublas)
 * @details Returns the norm() of the Eigen expression.
 * @tparam TDerived Derived type of the Eigen expression.
 * @param rM Vector or matrix.
 * @return Frobenius norm of rM.
 */
template<class TDerived>
inline auto norm_frobenius(const Eigen::MatrixBase<TDerived>& rM)
{
    return rM.norm();
}

/**
 * @brief Element-wise (Hadamard) product
 * @details Returns the lazy Eigen expression rX.cwiseProduct(rY).
 * @tparam TDerived1 Derived type of the first Eigen expression.
 * @tparam TDerived2 Derived type of the second Eigen expression.
 * @param rX First operand.
 * @param rY Second operand.
 * @return Lazy Eigen expression of the element-wise product of rX and rY.
 */
template<class TDerived1, class TDerived2>
inline auto element_prod(const Eigen::MatrixBase<TDerived1>& rX, const Eigen::MatrixBase<TDerived2>& rY)
{
    return rX.derived().cwiseProduct(rY.derived());
}

/**
 * @brief Element-wise division
 * @details Returns the lazy Eigen expression rX.cwiseQuotient(rY).
 * @tparam TDerived1 Derived type of the first Eigen expression.
 * @tparam TDerived2 Derived type of the second Eigen expression.
 * @param rX Dividend, the expression whose entries are divided.
 * @param rY Divisor, the expression whose entries divide those of rX.
 * @return Lazy Eigen expression of the element-wise quotient of rX by rY.
 */
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
 * @brief noalias() proxy: assignment without a protective temporary
 * @details Dense assignments go through Eigen's own NoAlias, with the ublas element-wise scalar conversion (the cast is
 * the identity for matching scalars), while sparse assignments use the plain Eigen operators.
 * @tparam TTargetHolder Target expression type, held by reference (lvalue targets) or by value (temporary proxies
 * such as row(M, i)).
 */
template<class TTargetHolder>
class NoAliasProxy
{
public:
    /// Type of the target expression. The holder type without its reference.
    using TargetType = std::remove_reference_t<TTargetHolder>;

    /// Scalar type of the target. Scalar type the dense expressions are cast to before being assigned.
    using Scalar = typename TargetType::Scalar;

    /**
     * @brief Constructs the proxy from an lvalue target
     * @details The target is bound to the holder as is, so it is referenced when TTargetHolder is a reference type.
     * @param rTarget Target of the assignments, referenced when TTargetHolder is a reference type.
     */
    explicit NoAliasProxy(TargetType& rTarget) : mTarget(rTarget) {}

    /**
     * @brief Constructs the proxy from a temporary target
     * @details The target is moved into the holder, which is then a value type (temporary proxies such as row(M, i)).
     * @param rTarget Temporary target, moved into the holder.
     */
    explicit NoAliasProxy(TargetType&& rTarget) : mTarget(std::move(rTarget)) {}

    /**
     * @brief Assigns a dense expression to the target
     * @details The expression is cast to the scalar type of the target (the identity for matching scalars) and assigned
     * through Eigen's noalias(), so no protective temporary is created. Returns the target.
     * @tparam TDerived Derived type of the dense Eigen expression.
     * @param rExpression Dense expression assigned to the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() = rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    /**
     * @brief Adds a dense expression to the target
     * @details The expression is cast to the scalar type of the target (the identity for matching scalars) and added
     * through Eigen's noalias(), so no protective temporary is created. Returns the target.
     * @tparam TDerived Derived type of the dense Eigen expression.
     * @param rExpression Dense expression added to the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() += rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    /**
     * @brief Subtracts a dense expression from the target
     * @details The expression is cast to the scalar type of the target (the identity for matching scalars) and
     * subtracted through Eigen's noalias(), so no protective temporary is created. Returns the target.
     * @tparam TDerived Derived type of the dense Eigen expression.
     * @param rExpression Dense expression subtracted from the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mTarget.noalias() -= rExpression.derived().template cast<Scalar>();
        return mTarget;
    }

    /**
     * @brief Assigns a sparse expression to the target
     * @details It uses the plain Eigen assignment, with neither a scalar cast nor noalias(). Returns the target.
     * @tparam TDerived Derived type of the sparse Eigen expression.
     * @param rExpression Sparse expression assigned to the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget = rExpression.derived();
        return mTarget;
    }

    /**
     * @brief Adds a sparse expression to the target
     * @details It uses the plain Eigen addition assignment, with neither a scalar cast nor noalias(). Returns the
     * target.
     * @tparam TDerived Derived type of the sparse Eigen expression.
     * @param rExpression Sparse expression added to the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator+=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget += rExpression.derived();
        return mTarget;
    }

    /**
     * @brief Subtracts a sparse expression from the target
     * @details It uses the plain Eigen subtraction assignment, with neither a scalar cast nor noalias(). Returns the
     * target.
     * @tparam TDerived Derived type of the sparse Eigen expression.
     * @param rExpression Sparse expression subtracted from the target.
     * @return Reference to the target.
     */
    template<class TDerived>
    TargetType& operator-=(const Eigen::SparseMatrixBase<TDerived>& rExpression)
    {
        mTarget -= rExpression.derived();
        return mTarget;
    }

private:
    /**
     * @brief Target of the assignments
     * @details Held by reference or by value, as selected by TTargetHolder.
     */
    TTargetHolder mTarget;
};
} // namespace Internals

/**
 * @brief noalias() on a dense lvalue target (a container or a stored proxy)
 * @details Returns an Internals::NoAliasProxy referencing the target. The assignments made through it skip the
 * protective temporary.
 * @tparam TDerived Derived type of the dense target.
 * @param rTarget Dense lvalue target of the assignments.
 * @return Internals::NoAliasProxy referencing the target.
 */
template<class TDerived>
inline auto noalias(Eigen::MatrixBase<TDerived>& rTarget)
{
    return Internals::NoAliasProxy<TDerived&>(rTarget.derived());
}

/**
 * @brief noalias() on a temporary dense proxy, as in noalias(row(M, i)) = expression
 * @details Returns an Internals::NoAliasProxy that holds the temporary by value (it is moved into the proxy).
 * @tparam TDerived Derived type of the dense target.
 * @param rTarget Temporary dense proxy, moved into the returned proxy.
 * @return Internals::NoAliasProxy holding the temporary by value.
 */
template<class TDerived>
inline auto noalias(Eigen::MatrixBase<TDerived>&& rTarget)
{
    return Internals::NoAliasProxy<TDerived>(std::move(rTarget.derived()));
}

/**
 * @brief noalias() on a sparse lvalue target
 * @details Returns an Internals::NoAliasProxy referencing the target. The sparse assignments of the proxy use the plain
 * Eigen operators, with no scalar cast.
 * @tparam TDerived Derived type of the sparse target.
 * @param rTarget Sparse lvalue target of the assignments.
 * @return Internals::NoAliasProxy referencing the target.
 */
template<class TDerived>
inline auto noalias(Eigen::SparseMatrixBase<TDerived>& rTarget)
{
    return Internals::NoAliasProxy<TDerived&>(rTarget.derived());
}

///@}

} // namespace Kratos

namespace Eigen
{
/**
 * @brief ADL hook for the member overload noalias(target) injected into Eigen::MatrixBase
 * @details The member overload is injected by includes/eigen_matrixbase_plugin.h. Inside a class deriving from an Eigen
 * type, an unqualified noalias(x) finds Eigen's member noalias() and never reaches Kratos::noalias. The plugin adds a
 * one-argument member overload that calls this hook unqualified; it lives in namespace Eigen because every Eigen-based
 * target has Eigen as an associated namespace, so ADL always finds it. The hook forwards its argument to
 * Kratos::noalias.
 * @tparam TTarget Type of the target, deduced as a forwarding reference.
 * @param rTarget Target forwarded to Kratos::noalias().
 * @return Proxy returned by Kratos::noalias() for the target.
 */
template<class TTarget>
inline auto KratosNoAliasHook(TTarget&& rTarget)
{
    return Kratos::noalias(std::forward<TTarget>(rTarget));
}
} // namespace Eigen

namespace Kratos
{

///@name LU factorization (boost::numeric::ublas::lu_factorize / lu_substitute)
///@{

/**
 * @brief Permutation (transposition sequence) produced by lu_factorize, as boost::numeric::ublas::permutation_matrix
 * @details The entry i holds the row exchanged with the row i at the step i of the factorization, and i itself when no
 * exchange was needed. ApplyTo() replays the sequence on the rows of a vector or matrix.
 * @tparam T Integer type of the stored indices, std::size_t by default.
 */
template<class T = std::size_t>
class permutation_matrix : public EigenVector<T>
{
public:
    /// Vector type this class derives from. Dense Eigen-backed vector holding the indices of the transposition
    /// sequence.
    using BaseType = EigenVector<T>;

    /// Type of the sizes of the permutation. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Constructs the identity permutation of the given size
     * @details The entry i is initialized to i, that is no row is exchanged.
     * @param Size Number of entries of the permutation.
     */
    explicit permutation_matrix(const std::size_t Size) : BaseType(Size)
    {
        for (std::size_t i = 0; i < Size; ++i) {
            (*this)[i] = static_cast<T>(i);
        }
    }

    /**
     * @brief Applies the transposition sequence to the rows of a vector or matrix (ublas swap_rows)
     * @details For each position i, in increasing order, the rows i and (*this)[i] of rX are swapped when they differ.
     * @tparam TDerived Derived type of the dense Eigen expression.
     * @param rX Vector or matrix whose rows are exchanged in place.
     */
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
 * @brief In-place LU factorization with partial pivoting, as boost::numeric::ublas::lu_factorize(m, pm)
 * @details A direct port of the uBLAS algorithm, so the packed L and U layout, the
 * transposition sequence written into rPM (used e.g. by the determinant sign
 * in MathUtils) and the return value (0, or i + 1 for the first zero pivot)
 * are identical to the uBLAS ones.
 * @tparam TDerived Derived type of the dense matrix.
 * @tparam TPermutation Type of the permutation, written through operator()(std::size_t), for example
 * permutation_matrix.
 * @param rM Matrix to factorize, overwritten in place with the packed L and U factors.
 * @param rPM Permutation receiving the transposition sequence of the pivoting: the entry i is set to the row exchanged
 * with the row i when a pivot swap happens; entries of the steps without a swap are left unchanged.
 * @return 0 if no zero pivot was found, otherwise i + 1 for the first zero pivot, found at the step i.
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

/**
 * @brief Forward/back substitution with the factors of lu_factorize, as boost::numeric::ublas::lu_substitute
 * @details The call mirrors lu_substitute(m, pm, x). rX (a vector or a matrix of right-hand sides) is overwritten with
 * the solution: the row exchanges of the permutation are applied first, then the unit lower and the upper triangular
 * factors stored in rM are solved in place.
 * @tparam TDerived Derived type of the matrix holding the factors.
 * @tparam TPermutation Type of the permutation, read through size() and operator()(std::size_t), for example
 * permutation_matrix.
 * @tparam TDerivedRhs Derived type of the right-hand side.
 * @param rM Packed L and U factors, as produced by lu_factorize().
 * @param rPM Transposition sequence produced by lu_factorize().
 * @param rX Right-hand side, a vector or a matrix of right-hand sides, overwritten with the solution.
 */
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

/// Index array of an indirect view. Any indexable container works (ublas indirect_array<A>): the alias is the array
/// type itself.
template<class TArray>
using indirect_array = TArray;

namespace Internals
{
/**
 * @brief Nullary functor reading a matrix through two index arrays
 * @details For the indices (I, J) of the view it returns the entry of the matrix at the row mpRowIndices[I] and the
 * column mpColumnIndices[J]. It is the operation of the Eigen nullary expression of matrix_indirect, and it refers to
 * the matrix and to the index arrays through pointers.
 * @tparam TMatrix Type of the matrix that is read.
 * @tparam TIndexArray Type of the row and column index arrays, indexable by std::size_t.
 */
template<class TMatrix, class TIndexArray>
struct IndirectAccessOp
{
    /**
     * @brief Matrix that is read
     * @details Referenced through a pointer, not owned.
     */
    const TMatrix* mpMatrix;

    /**
     * @brief Row index array
     * @details The entry I of the array is the row of the matrix read for the row I of the view.
     */
    const TIndexArray* mpRowIndices;

    /**
     * @brief Column index array
     * @details The entry J of the array is the column of the matrix read for the column J of the view.
     */
    const TIndexArray* mpColumnIndices;

    /**
     * @brief Entry (I, J) of the indirect view
     * @details Returns the entry of the matrix at the row mpRowIndices[I] and the column mpColumnIndices[J].
     * @param I Row of the view.
     * @param J Column of the view.
     * @return Entry of the matrix read through the row and column index arrays.
     */
    typename TMatrix::value_type operator()(const Eigen::Index I, const Eigen::Index J) const
    {
        return (*mpMatrix)((*mpRowIndices)[static_cast<std::size_t>(I)], (*mpColumnIndices)[static_cast<std::size_t>(J)]);
    }
};
} // namespace Internals

/**
 * @brief Read-only indirect (gather) view of a matrix through two index arrays
 * @details As boost::numeric::ublas::matrix_indirect. It is a genuine Eigen expression, so it can be evaluated into any
 * dense type. The matrix and the index arrays are referenced, not copied, so they must outlive the view.
 * @tparam TMatrix Type of the matrix that is read, possibly const-qualified.
 * @tparam TIndexArray Type of the row and column index arrays, any container with size() indexable by std::size_t.
 */
template<class TMatrix, class TIndexArray>
class matrix_indirect : public Eigen::CwiseNullaryOp<Internals::IndirectAccessOp<std::remove_const_t<TMatrix>, TIndexArray>, Internals::DynamicRowMajorMatrix<typename std::remove_const_t<TMatrix>::value_type>>
{
public:
    /// Operation type of the nullary expression. Internals::IndirectAccessOp over the matrix type without its const
    /// qualifier.
    using OpType = Internals::IndirectAccessOp<std::remove_const_t<TMatrix>, TIndexArray>;

    /// Eigen nullary expression this view derives from. A dynamic row-major matrix expression whose entries are
    /// computed by the operation type.
    using BaseType = Eigen::CwiseNullaryOp<OpType, Internals::DynamicRowMajorMatrix<typename std::remove_const_t<TMatrix>::value_type>>;

    /// Scalar type of the entries of the view. Taken from the matrix type without its const qualifier.
    using value_type = typename std::remove_const_t<TMatrix>::value_type;

    /// Type of the sizes of the view. Alias of std::size_t, as in the uBLAS surface.
    using size_type = std::size_t;

    /**
     * @brief Constructs the view of a matrix through a row and a column index array
     * @details The view has as many rows as rRowIndices and as many columns as rColumnIndices, and its entry (I, J) is
     * rMatrix(rRowIndices[I], rColumnIndices[J]).
     * @param rMatrix Matrix that is read, referenced and not copied.
     * @param rRowIndices Row index array, giving the row of the matrix read for each row of the view.
     * @param rColumnIndices Column index array, giving the column of the matrix read for each column of the view.
     */
    matrix_indirect(const TMatrix& rMatrix, const TIndexArray& rRowIndices, const TIndexArray& rColumnIndices)
        : BaseType(static_cast<Eigen::Index>(rRowIndices.size()), static_cast<Eigen::Index>(rColumnIndices.size()), OpType{&rMatrix, &rRowIndices, &rColumnIndices}) {}

    /**
     * @brief Number of rows of the view
     * @details Returns the number of rows of the Eigen base expression as an unsigned integer, as in the uBLAS surface.
     * @return Number of rows of the view.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns of the view
     * @details Returns the number of columns of the Eigen base expression as an unsigned integer, as in the uBLAS
     * surface.
     * @return Number of columns of the view.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

///@}

} // namespace Kratos
