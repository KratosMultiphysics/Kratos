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
#include <cstddef>
#include <istream>
#include <iterator>
#include <ostream>
#include <type_traits>
#include <utility>
#include <memory>

// External includes
#include <Eigen/Core>

// Project includes
#include "includes/define.h"
#include "includes/storage_view.h"

namespace Kratos
{

///@name Kratos Classes
///@{

namespace Internals
{

/**
 * @brief Compile-time flag for a column-vector-shaped Eigen expression (N x 1)
 * @details True when the expression has one column and not one row at compile time, the shape uBLAS treats as a
 * vector type.
 * @tparam TDerived Eigen expression type whose shape is tested.
 */
template<class TDerived>
inline constexpr bool IsColumnShaped = TDerived::ColsAtCompileTime == 1 && TDerived::RowsAtCompileTime != 1;

/**
 * @brief Compile-time flag for a vector-shaped Eigen expression of either orientation
 * @details True when the expression has one column or one row at compile time.
 * @tparam TDerived Eigen expression type whose shape is tested.
 */
template<class TDerived>
inline constexpr bool IsVectorShaped = TDerived::ColsAtCompileTime == 1 || TDerived::RowsAtCompileTime == 1;

/**
 * @brief Static-shape compatibility of an Eigen expression with a TSize1 x TSize2 bounded matrix
 * @details True when each compile-time extent of the expression is either dynamic or equal to the bounded one.
 * Vector shapes of either orientation are compatible with a bounded row/column matrix of the same length (Eigen
 * transposes vectors on assignment), as uBLAS assigns a transposed vector into a 1 x N matrix.
 * @tparam TDerived Eigen expression type whose shape is tested.
 * @tparam TSize1 Static number of rows of the bounded matrix.
 * @tparam TSize2 Static number of columns of the bounded matrix.
 */
template<class TDerived, std::size_t TSize1, std::size_t TSize2>
inline constexpr bool IsBoundedCompatible =
    ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
     (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2))) ||
    (TSize1 == 1 && TDerived::ColsAtCompileTime == 1 && (TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize2))) ||
    (TSize2 == 1 && TDerived::RowsAtCompileTime == 1 && (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize1)));

/**
 * @brief Storage options for the Eigen-backed bounded matrix
 * @details Row-major to match the uBLAS bounded_matrix memory layout, except where Eigen mandates the orientation
 * (column vectors must be column-major, row vectors row-major): the result is Eigen::ColMajor when Size2 is 1 and
 * Size1 is not, Eigen::RowMajor otherwise, always combined with Eigen::DontAlign. DontAlign keeps the plain
 * T[N1*N2] layout without over-alignment. Size1 and Size2 are the static extents of the bounded matrix.
 * @param Size1 Static number of rows of the bounded matrix.
 * @param Size2 Static number of columns of the bounded matrix.
 * @return Eigen storage options: Eigen::ColMajor for a Size1 x 1 shape with Size1 != 1, Eigen::RowMajor otherwise,
 * combined with Eigen::DontAlign.
 */
constexpr int EigenBoundedMatrixOptions(const std::size_t Size1, const std::size_t Size2)
{
    return ((Size2 == 1 && Size1 != 1) ? Eigen::ColMajor : Eigen::RowMajor) | Eigen::DontAlign;
}

/// Dynamic Eigen matrix of the given storage order (Eigen::RowMajor or Eigen::ColMajor). The plain object behind the
/// Kratos dense matrix types. Row-major by default, matching uBLAS's own default (row_major) and the storage order of
/// the uBLAS dense matrices, so raw-buffer interoperability is preserved; an explicit Eigen::ColMajor is used where a
/// routine written against a foreign (Fortran/LAPACK-style) API needs column-major storage (see Kratos::column_major in
/// eigen_interface.h). T is the scalar type and TLayout the storage order.
template<class T, int TLayout = Eigen::RowMajor>
using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, TLayout>;

/// Row-major dynamic Eigen matrix. Kept for code written against the previous (row-major only) name. T is the scalar
/// type.
template<class T>
using DynamicRowMajorMatrix = DynamicMatrix<T, Eigen::RowMajor>;

/// Dynamic Eigen column vector. The plain object behind the Kratos dense vectors. T is the scalar type.
template<class T>
using DynamicColumnVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

} // namespace Internals

/**
 * @brief Dense, dynamically sized Eigen matrix with the uBLAS member surface
 * @details Row-major by default, matching the storage order of the uBLAS dense matrices (so raw-buffer
 * interoperability is preserved); the second template parameter lets a handful of call sites request column-major
 * storage instead (e.g. a local buffer handed to a Fortran/LAPACK-style routine, matching uBLAS's own
 * matrix<T, column_major>). The added members (size1/size2, the preserving resize, clear, data(), the assignment
 * protocol) mirror boost::numeric::ublas::matrix so generic Kratos code compiles unchanged. Eigen only ever sees the
 * base class, so shadowing its members here is safe.
 * @tparam TDataType Scalar type stored in the matrix (e.g. double).
 * @tparam TLayout Eigen storage order: Eigen::RowMajor (the default) or Eigen::ColMajor.
 */
template<class TDataType, int TLayout = Eigen::RowMajor>
class EigenMatrix : public Internals::DynamicMatrix<TDataType, TLayout>
{
public:
    ///@name Type Definitions
    ///@{

    /// Plain Eigen matrix type this class derives from. Dynamic-extent Eigen matrix of TDataType with the storage order
    /// TLayout.
    using BaseType = Internals::DynamicMatrix<TDataType, TLayout>;

    /// Scalar type of the stored entries. Alias of TDataType, as in the uBLAS matrix.
    using value_type = TDataType;

    /// Unsigned type used for sizes and indices. Alias of std::size_t, as in the uBLAS matrix.
    using size_type = std::size_t;

    /// Signed type used for differences between positions. Alias of std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;

    /// Mutable reference to a stored entry. Alias of TDataType&.
    using reference = TDataType&;

    /// Read-only reference to a stored entry. Alias of const TDataType&.
    using const_reference = const TDataType&;

    /// Pointer to a stored entry. Alias of TDataType*.
    using pointer = TDataType*;

    /// Iterator over the entries in storage order. Plain pointer, since the storage is contiguous (see begin() and
    /// end()).
    using iterator = TDataType*;

    /// Read-only iterator over the entries in storage order. Plain pointer to const, since the storage is contiguous
    /// (see begin() and end()).
    using const_iterator = const TDataType*;

    /// Type of the contiguous storage array exposed by data(). Internals::StorageView of the entries, a raw pointer
    /// plus the length that also converts implicitly to the raw pointer. It plays the role of the uBLAS storage array.
    using array_type = Internals::StorageView<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates an empty (0 x 0) matrix.
     */
    EigenMatrix() = default;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and every entry of rOther.
     * @param rOther Matrix to copy.
     */
    EigenMatrix(const EigenMatrix& rOther) = default;

    /**
     * @brief Move constructor
     * @details Defaulted: takes over the storage of rOther.
     * @param rOther Matrix whose storage is taken over.
     */
    EigenMatrix(EigenMatrix&& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Defaulted: copies the extents and every entry of rOther.
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    EigenMatrix& operator=(const EigenMatrix& rOther) = default;

    /**
     * @brief Move assignment operator
     * @details Defaulted: takes over the storage of rOther.
     * @param rOther Matrix whose storage is taken over.
     * @return Reference to this matrix.
     */
    EigenMatrix& operator=(EigenMatrix&& rOther) = default;

    /**
     * @brief Constructor with the extents
     * @details Allocates an uninitialized Size1 x Size2 matrix.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    EigenMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /**
     * @brief Fill constructor
     * @details uBLAS-style (rows, cols, value) constructor: allocates a Size1 x Size2 matrix and sets every entry to
     * rValue.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     * @param rValue Value assigned to every entry.
     */
    EigenMatrix(const std::size_t Size1, const std::size_t Size2, const TDataType& rValue) : BaseType(Size1, Size2)
    {
        std::fill(begin(), end(), rValue);
    }

    /**
     * @brief Constructor from any matrix-shaped Eigen expression
     * @details The scalar cast makes the cross-precision conversions of ublas work; it is the identity for matching
     * scalars. Column-vector-shaped (N x 1) expressions are excluded, as ublas does not construct a matrix from a
     * vector expression either.
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression whose values, cast to TDataType, initialize the matrix.
     */
    template<class TDerived>
    requires (!Internals::IsColumnShaped<TDerived>)
    EigenMatrix(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment from any matrix-shaped Eigen expression
     * @details The assignment evaluates through a temporary (.eval()), reproducing the alias-safe uBLAS operator=
     * semantics (e.g. M = trans(M)); the temporary-free fast path remains noalias(M) = expr or assign(), as in uBLAS.
     * Column-vector-shaped (N x 1) expressions are excluded, as in the constructor.
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to assign to this matrix, cast to TDataType.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    requires (!Internals::IsColumnShaped<TDerived>)
    EigenMatrix& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>().eval());
        return *this;
    }

    // Element access. ublas::matrix narrows whatever index it is given to its
    // size_type, so m(i, j) reads a coefficient even when i/j are the doubles
    // that come out of a Vector; Eigen instead reads a pair of floating-point
    // arguments as a fancy IndexedView (its "A(indices, indices)" indexing).
    // The first overload restores the ublas reading for any pair of arithmetic
    // indices, the other two forward everything else (index arrays, Eigen::all,
    // seq, the single-argument linear access) to Eigen unchanged. They are
    // written as one complete set rather than added next to a
    // "using BaseType::operator();": mixing the two makes the calls with
    // converted indices ambiguous.

    /**
     * @brief uBLAS-style coefficient access for any pair of arithmetic indices
     * @details The indices are converted to Eigen::Index, so m(i, j) reads a coefficient even when i and j are
     * floating-point values. Returns a reference to the entry.
     * @tparam TIndex1 Arithmetic type of the row index.
     * @tparam TIndex2 Arithmetic type of the column index.
     * @param I Row index.
     * @param J Column index.
     * @return Reference to the entry at row I and column J.
     */
    template<class TIndex1, class TIndex2>
    requires (std::is_arithmetic_v<TIndex1> && std::is_arithmetic_v<TIndex2>)
    reference operator()(const TIndex1 I, const TIndex2 J)
    {
        return BaseType::operator()(static_cast<Eigen::Index>(I), static_cast<Eigen::Index>(J));
    }

    /**
     * @brief uBLAS-style coefficient access for any pair of arithmetic indices (const version)
     * @details Same as the mutable overload, but returns a const reference to the entry.
     * @tparam TIndex1 Arithmetic type of the row index.
     * @tparam TIndex2 Arithmetic type of the column index.
     * @param I Row index.
     * @param J Column index.
     * @return Const reference to the entry at row I and column J.
     */
    template<class TIndex1, class TIndex2>
    requires (std::is_arithmetic_v<TIndex1> && std::is_arithmetic_v<TIndex2>)
    const_reference operator()(const TIndex1 I, const TIndex2 J) const
    {
        return BaseType::operator()(static_cast<Eigen::Index>(I), static_cast<Eigen::Index>(J));
    }

    /**
     * @brief Two-argument subscript other than a pair of arithmetic indices
     * @details Everything else Eigen accepts as a two-argument subscript (index arrays, Eigen::all, seq, ...),
     * forwarded unchanged to Eigen.
     * @tparam TArg1 Type of the row subscript, anything Eigen accepts as a row selector.
     * @tparam TArg2 Type of the column subscript, anything Eigen accepts as a column selector.
     * @param rArgument1 Row subscript, forwarded unchanged to Eigen.
     * @param rArgument2 Column subscript, forwarded unchanged to Eigen.
     * @return Eigen's result for this subscript: a coefficient, a block or an indexed view, depending on the arguments.
     */
    template<class TArg1, class TArg2>
    requires (!std::is_arithmetic_v<std::decay_t<TArg1>> || !std::is_arithmetic_v<std::decay_t<TArg2>>)
    decltype(auto) operator()(TArg1&& rArgument1, TArg2&& rArgument2)
    {
        return BaseType::operator()(std::forward<TArg1>(rArgument1), std::forward<TArg2>(rArgument2));
    }

    /**
     * @brief Two-argument subscript other than a pair of arithmetic indices (const version)
     * @details Forwarded unchanged to Eigen.
     * @tparam TArg1 Type of the row subscript, anything Eigen accepts as a row selector.
     * @tparam TArg2 Type of the column subscript, anything Eigen accepts as a column selector.
     * @param rArgument1 Row subscript, forwarded unchanged to Eigen.
     * @param rArgument2 Column subscript, forwarded unchanged to Eigen.
     * @return Eigen's read-only result for this subscript: a coefficient, a block or an indexed view.
     */
    template<class TArg1, class TArg2>
    requires (!std::is_arithmetic_v<std::decay_t<TArg1>> || !std::is_arithmetic_v<std::decay_t<TArg2>>)
    decltype(auto) operator()(TArg1&& rArgument1, TArg2&& rArgument2) const
    {
        return BaseType::operator()(std::forward<TArg1>(rArgument1), std::forward<TArg2>(rArgument2));
    }

    /**
     * @brief Single-argument subscript (Eigen's linear access)
     * @details Forwarded unchanged to Eigen.
     * @tparam TArg Type of the subscript argument, deduced from the call.
     * @param rArgument Linear index of the entry, forwarded unchanged to Eigen.
     * @return Eigen's result for this subscript: a reference to the entry at the linear index.
     */
    template<class TArg>
    decltype(auto) operator()(TArg&& rArgument)
    {
        return BaseType::operator()(std::forward<TArg>(rArgument));
    }

    /**
     * @brief Single-argument subscript (const version)
     * @details Eigen's linear access, forwarded unchanged to Eigen.
     * @tparam TArg Type of the subscript argument, deduced from the call.
     * @param rArgument Linear index of the entry, forwarded unchanged to Eigen.
     * @return Eigen's result for this subscript: the entry at the linear index, read-only.
     */
    template<class TArg>
    decltype(auto) operator()(TArg&& rArgument) const
    {
        return BaseType::operator()(std::forward<TArg>(rArgument));
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief uBLAS assignment protocol: assignment without a protective temporary
     * @details Resizes the target when needed (a superset of the ublas contract). The expression is cast to TDataType
     * and assigned through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to assign to this matrix, cast to TDataType.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenMatrix& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place addition without a protective temporary
     * @details The expression is cast to TDataType and added to this matrix through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to add to this matrix, cast to TDataType.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenMatrix& plus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() += rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place subtraction without a protective temporary
     * @details The expression is cast to TDataType and subtracted from this matrix through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to subtract from this matrix, cast to TDataType.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenMatrix& minus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() -= rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS-style resize
     * @details As in ublas::matrix the default preserves the values: when Preserve is true the entries common to the
     * old and the new extents are kept, using Eigen's conservativeResize or, for move-only scalars (e.g. a container
     * of sparse graphs), a move into freshly allocated storage because conservativeResize copy-assigns. When Preserve
     * is false the matrix is resized without keeping its values.
     * @param NewSize1 New number of rows.
     * @param NewSize2 New number of columns.
     * @param Preserve If true (the default) the entries common to the old and the new extents are kept.
     */
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) {
            if constexpr (std::is_copy_assignable_v<TDataType>) {
                this->conservativeResize(NewSize1, NewSize2);
            } else {
                // Move-only scalars (e.g. a container of sparse graphs): Eigen's
                // conservativeResize copy-assigns, so keep the entries by moving
                // them into freshly allocated storage instead.
                BaseType new_storage(NewSize1, NewSize2);
                const std::size_t common_rows = std::min(NewSize1, size1());
                const std::size_t common_cols = std::min(NewSize2, size2());
                for (std::size_t i = 0; i < common_rows; ++i) {
                    for (std::size_t j = 0; j < common_cols; ++j) {
                        new_storage(i, j) = std::move((*this)(i, j));
                    }
                }
                BaseType::swap(new_storage);
            }
        } else {
            BaseType::resize(NewSize1, NewSize2);
        }
    }

    /**
     * @brief uBLAS-style clear
     * @details Value-initializes every entry; the size is kept.
     */
    void clear()
    {
        std::fill(begin(), end(), TDataType());
    }

    /**
     * @brief Swaps the contents with another matrix
     * @details Delegates to Eigen's swap, which for dynamic matrices of the same type only exchanges the data
     * pointers.
     * @param rOther Matrix whose contents are exchanged with this one.
     */
    void swap(EigenMatrix& rOther)
    {
        BaseType::swap(rOther);
    }

    /**
     * @brief Swaps the contents of two matrices
     * @details Non-member overload found by argument-dependent lookup; it forwards to the swap member.
     * @param rFirst Matrix whose contents are exchanged with rSecond.
     * @param rSecond Matrix whose contents are exchanged with rFirst.
     */
    friend void swap(EigenMatrix& rFirst, EigenMatrix& rSecond)
    {
        rFirst.swap(rSecond);
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Number of rows
     * @details Returns the row count as an unsigned value (Eigen's rows() is signed).
     * @return Number of rows.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns
     * @details Returns the column count as an unsigned value (Eigen's cols() is signed).
     * @return Number of columns.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /**
     * @brief Contiguous storage view of the entries (in TLayout order)
     * @details Provides the uBLAS data() surface and an implicit conversion to the raw pointer.
     * @return Storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<TDataType> data() { return {BaseType::data(), static_cast<std::size_t>(BaseType::size())}; }

    /**
     * @brief Contiguous read-only storage view of the entries (in TLayout order)
     * @details Const version of data(): the view refers to const entries.
     * @return Read-only storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<const TDataType> data() const { return {BaseType::data(), static_cast<std::size_t>(BaseType::size())}; }

    /**
     * @brief Iterator to the first entry
     * @details Pointer to the contiguous storage, traversed in TLayout order (row by row for the default row-major
     * layout).
     * @return Pointer to the first entry.
     */
    iterator begin() { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry
     * @details Const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator begin() const { return BaseType::data(); }

    /**
     * @brief Iterator past the last entry
     * @details Pointer one past the end of the contiguous storage.
     * @return Pointer one past the last entry.
     */
    iterator end() { return BaseType::data() + BaseType::size(); }

    /**
     * @brief Read-only iterator past the last entry
     * @details Const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator end() const { return BaseType::data() + BaseType::size(); }

    ///@}
    ///@name Row iterators
    ///@{
    // Row-major traversal with the boost::numeric::ublas::matrix iterator
    // concept, so the ublas idiom compiles unchanged:
    //     for (auto i1 = m.begin1(); i1 != m.end1(); ++i1)
    //         for (auto i2 = i1.begin(); i2 != i1.end(); ++i2) ...
    //     std::copy(v.begin(), v.end(), (m.begin1() + i).begin());
    // A row of the (row-major) storage is contiguous, so the entry iterators
    // (iterator2) are plain pointers.

    /// Iterator over the entries of one row. Plain pointer, since a row of the row-major storage is contiguous.
    using iterator2 = TDataType*;

    /// Read-only iterator over the entries of one row. Plain pointer to const, since a row of the row-major storage is
    /// contiguous.
    using const_iterator2 = const TDataType*;

    /**
     * @brief Iterator over the rows of the matrix (ublas iterator1)
     * @details begin()/end() yield the entries of the current row.
     * @tparam TIsConst Selects the const flavour: when true the matrix is held through a pointer to const and the
     * entry iterators are read-only, otherwise they grant write access.
     */
    template<bool TIsConst>
    class RowIterator
    {
    public:
        /// Pointer type to the traversed matrix. Pointer to const EigenMatrix when TIsConst is true, pointer to
        /// EigenMatrix otherwise.
        using MatrixPointerType = std::conditional_t<TIsConst, const EigenMatrix*, EigenMatrix*>;

        /// Iterator type over the entries of a row. Pointer to const TDataType when TIsConst is true, pointer to
        /// TDataType otherwise.
        using EntryIteratorType = std::conditional_t<TIsConst, const TDataType*, TDataType*>;

        /**
         * @brief Constructor from a matrix and a row
         * @details Creates an iterator pointing at the row Row of the matrix pointed to by pMatrix.
         * @param pMatrix Pointer to the matrix to traverse (pointer to const for the const iterator).
         * @param Row Index of the row pointed at.
         */
        RowIterator(MatrixPointerType pMatrix, const std::size_t Row) : mpMatrix(pMatrix), mRow(Row) {}

        /**
         * @brief Conversion from the mutable to the const row iterator
         * @details Same conversion as in ublas: only the mutable-to-const direction is available, so it requires
         * TIsConst to be true and the source iterator to be the mutable one.
         * @tparam TOtherIsConst Constness of the source iterator; the constraint requires it to be false (mutable).
         * @param rOther Mutable row iterator to convert from.
         */
        template<bool TOtherIsConst>
        requires (TIsConst && !TOtherIsConst)
        RowIterator(const RowIterator<TOtherIsConst>& rOther) : mpMatrix(rOther.mpMatrix), mRow(rOther.mRow) {}

        /**
         * @brief Index of the row currently pointed at
         * @details Returns the row index held by the iterator.
         * @return Index of the row currently pointed at.
         */
        std::size_t index1() const { return mRow; }

        /**
         * @brief First entry of this row
         * @details Pointer to the first entry of the current row, computed on the contiguous row-major storage as the
         * row index times the number of columns.
         * @return Pointer to the first entry of the current row.
         */
        EntryIteratorType begin() const { return mpMatrix->BaseType::data() + mRow * mpMatrix->size2(); }

        /**
         * @brief Past-the-last entry of this row
         * @details Pointer one past the last entry of the current row, that is begin() plus the number of columns.
         * @return Pointer one past the last entry of the current row.
         */
        EntryIteratorType end() const { return begin() + mpMatrix->size2(); }

        /**
         * @brief Pre-increment operator
         * @details Moves to the next row and returns this iterator.
         * @return Reference to this iterator, moved to the next row.
         */
        RowIterator& operator++() { ++mRow; return *this; }

        /**
         * @brief Post-increment operator
         * @details Moves to the next row and returns a copy of the iterator as it was before the move.
         * @return Copy of the iterator as it was before the move.
         */
        RowIterator operator++(int) { RowIterator copy(*this); ++mRow; return copy; }

        /**
         * @brief Pre-decrement operator
         * @details Moves to the previous row and returns this iterator.
         * @return Reference to this iterator, moved to the previous row.
         */
        RowIterator& operator--() { --mRow; return *this; }

        /**
         * @brief Compound addition operator
         * @details Advances this iterator by Offset rows and returns it.
         * @param Offset Signed number of rows to advance by.
         * @return Reference to this iterator, advanced by Offset rows.
         */
        RowIterator& operator+=(const std::ptrdiff_t Offset) { mRow += Offset; return *this; }

        /**
         * @brief Addition operator
         * @details Returns a new iterator advanced by Offset rows; this iterator is left unchanged.
         * @param Offset Signed number of rows to advance by.
         * @return New iterator advanced by Offset rows.
         */
        RowIterator operator+(const std::ptrdiff_t Offset) const { return RowIterator(mpMatrix, mRow + Offset); }

        /**
         * @brief Difference between two row iterators
         * @details Returns the signed number of rows from rOther to this iterator (the difference of the row indices).
         * @param rOther Row iterator to subtract from this one.
         * @return Signed number of rows from rOther to this iterator.
         */
        std::ptrdiff_t operator-(const RowIterator& rOther) const { return static_cast<std::ptrdiff_t>(mRow) - static_cast<std::ptrdiff_t>(rOther.mRow); }

        /**
         * @brief Equality comparison
         * @details Two iterators are equal when they point at the same row index; the matrix is not compared.
         * @param rOther Row iterator to compare with.
         * @return True if both iterators point at the same row index.
         */
        bool operator==(const RowIterator& rOther) const { return mRow == rOther.mRow; }

        /**
         * @brief Inequality comparison
         * @details Negation of operator==.
         * @param rOther Row iterator to compare with.
         * @return True if the iterators point at different row indices.
         */
        bool operator!=(const RowIterator& rOther) const { return !(*this == rOther); }

    private:
        template<bool TOtherIsConst> friend class RowIterator;

        /**
         * @brief Pointer to the matrix being traversed
         * @details Pointer to const for the const row iterator.
         */
        MatrixPointerType mpMatrix;

        /**
         * @brief Index of the row currently pointed at
         * @details Updated by the increment, decrement and compound addition operators.
         */
        std::size_t mRow;
    };

    /// Mutable iterator over the rows of the matrix. RowIterator whose entry iterators grant write access (ublas
    /// iterator1).
    using iterator1 = RowIterator<false>;

    /// Read-only iterator over the rows of the matrix. RowIterator whose entry iterators are read-only (ublas
    /// const_iterator1).
    using const_iterator1 = RowIterator<true>;

    /**
     * @brief First row
     * @details Returns a mutable row iterator pointing at the first row. Row iterators require the row-major layout,
     * which is enforced with a static_assert.
     * @return Mutable row iterator pointing at the first row.
     */
    iterator1 begin1()
    {
        static_assert(TLayout == Eigen::RowMajor, "Row iterators require the row-major layout.");
        return iterator1(this, 0);
    }

    /**
     * @brief Past-the-last row
     * @details Returns a mutable row iterator pointing past the last row (row index size1()). Row iterators require
     * the row-major layout, which is enforced with a static_assert.
     * @return Mutable row iterator pointing past the last row.
     */
    iterator1 end1()
    {
        static_assert(TLayout == Eigen::RowMajor, "Row iterators require the row-major layout.");
        return iterator1(this, size1());
    }

    /**
     * @brief First row (const version)
     * @details Returns a read-only row iterator pointing at the first row. Row iterators require the row-major
     * layout, which is enforced with a static_assert.
     * @return Read-only row iterator pointing at the first row.
     */
    const_iterator1 begin1() const
    {
        static_assert(TLayout == Eigen::RowMajor, "Row iterators require the row-major layout.");
        return const_iterator1(this, 0);
    }

    /**
     * @brief Past-the-last row (const version)
     * @details Returns a read-only row iterator pointing past the last row (row index size1()). Row iterators
     * require the row-major layout, which is enforced with a static_assert.
     * @return Read-only row iterator pointing past the last row.
     */
    const_iterator1 end1() const
    {
        static_assert(TLayout == Eigen::RowMajor, "Row iterators require the row-major layout.");
        return const_iterator1(this, size1());
    }

    ///@}
};

/**
 * @brief Dense, dynamically sized Eigen column vector with the uBLAS member surface
 * @details Adds the (size), (size, value) constructors, the preserving resize, the unsigned size(), data() and the
 * assignment protocol of boost::numeric::ublas::vector; operator[] and operator() are already provided by Eigen with
 * compatible semantics.
 * @tparam TDataType Scalar type stored in the vector (e.g. double).
 */
template<class TDataType>
class EigenVector : public Internals::DynamicColumnVector<TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Plain Eigen vector type this class derives from. Dynamic-extent Eigen column vector of TDataType.
    using BaseType = Internals::DynamicColumnVector<TDataType>;

    /// Scalar type of the stored entries. Alias of TDataType, as in the uBLAS vector.
    using value_type = TDataType;

    /// Unsigned type used for sizes and indices. Alias of std::size_t, as in the uBLAS vector.
    using size_type = std::size_t;

    /// Signed type used for differences between positions. Alias of std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;

    /// Mutable reference to a stored entry. Alias of TDataType&.
    using reference = TDataType&;

    /// Read-only reference to a stored entry. Alias of const TDataType&.
    using const_reference = const TDataType&;

    /// Pointer to a stored entry. Alias of TDataType*.
    using pointer = TDataType*;

    /// Iterator over the entries. Plain pointer, since the storage is contiguous (see begin() and end()).
    using iterator = TDataType*;

    /// Read-only iterator over the entries. Plain pointer to const, since the storage is contiguous (see begin() and
    /// end()).
    using const_iterator = const TDataType*;

    /// Reverse iterator over the entries. std::reverse_iterator of iterator (see rbegin() and rend()).
    using reverse_iterator = std::reverse_iterator<iterator>;

    /// Read-only reverse iterator over the entries. std::reverse_iterator of const_iterator (see rbegin() and rend()).
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    /// Type of the contiguous storage array exposed by data(). Internals::StorageView of the entries, a raw pointer
    /// plus the length that also converts implicitly to the raw pointer. It plays the role of the uBLAS storage array.
    using array_type = Internals::StorageView<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates an empty (0-length) vector.
     */
    EigenVector() = default;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the size and every entry of rOther.
     * @param rOther Vector to copy.
     */
    EigenVector(const EigenVector& rOther) = default;

    /**
     * @brief Move constructor
     * @details Defaulted: takes over the storage of rOther.
     * @param rOther Vector whose storage is taken over.
     */
    EigenVector(EigenVector&& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Defaulted: copies the size and every entry of rOther.
     * @param rOther Vector to copy.
     * @return Reference to this vector.
     */
    EigenVector& operator=(const EigenVector& rOther) = default;

    /**
     * @brief Move assignment operator
     * @details Defaulted: takes over the storage of rOther.
     * @param rOther Vector whose storage is taken over.
     * @return Reference to this vector.
     */
    EigenVector& operator=(EigenVector&& rOther) = default;

    /**
     * @brief Constructor with the size
     * @details Allocates an uninitialized vector of the given size.
     * @param Size Number of entries.
     */
    explicit EigenVector(const std::size_t Size) : BaseType(Size) {}

    /**
     * @brief Fill constructor
     * @details uBLAS-style (size, value) constructor: allocates a vector of the given size and sets every entry to
     * rValue.
     * @param Size Number of entries.
     * @param rValue Value assigned to every entry.
     */
    EigenVector(const std::size_t Size, const TDataType& rValue) : BaseType(Size)
    {
        std::fill(begin(), end(), rValue);
    }

    /**
     * @brief Constructor from any vector-shaped Eigen expression
     * @details The scalar cast makes the cross-precision conversions of ublas work; it is the identity for matching
     * scalars. Vector shapes of either orientation are accepted.
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression whose values, cast to TDataType, initialize the vector.
     */
    template<class TDerived>
    requires (Internals::IsVectorShaped<TDerived>)
    EigenVector(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment from any vector-shaped Eigen expression
     * @details Evaluated through a temporary (alias-safe, as the ublas operator=); noalias(v) = expr and assign()
     * remain the temporary-free path. Vector shapes of either orientation are accepted.
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to assign to this vector, cast to TDataType.
     * @return Reference to this vector.
     */
    template<class TDerived>
    requires (Internals::IsVectorShaped<TDerived>)
    EigenVector& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>().eval());
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief uBLAS assignment protocol: assignment without a protective temporary
     * @details The expression is cast to TDataType and assigned through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to assign to this vector, cast to TDataType.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenVector& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place addition without a protective temporary
     * @details The expression is cast to TDataType and added to this vector through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to add to this vector, cast to TDataType.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenVector& plus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() += rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place subtraction without a protective temporary
     * @details The expression is cast to TDataType and subtracted from this vector through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rOther Expression to subtract from this vector, cast to TDataType.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenVector& minus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() -= rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS-style resize
     * @details As in ublas::vector the default preserves the values: when Preserve is true the leading entries common
     * to the old and the new size are kept, using Eigen's conservativeResize or, for move-only scalars (e.g. a
     * container of sparse graphs), a move into freshly allocated storage because conservativeResize copy-assigns.
     * When Preserve is false the vector is resized without keeping its values.
     * @param NewSize New number of entries.
     * @param Preserve If true (the default) the leading entries common to the old and the new size are kept.
     */
    void resize(const std::size_t NewSize, const bool Preserve = true)
    {
        if (Preserve) {
            if constexpr (std::is_copy_assignable_v<TDataType>) {
                this->conservativeResize(NewSize);
            } else {
                // Move-only scalars (e.g. a container of sparse graphs): Eigen's
                // conservativeResize copy-assigns, so keep the entries by moving
                // them into freshly allocated storage instead.
                BaseType new_storage(NewSize);
                const std::size_t common_size = std::min(NewSize, size());
                for (std::size_t i = 0; i < common_size; ++i) {
                    new_storage[i] = std::move((*this)[i]);
                }
                BaseType::swap(new_storage);
            }
        } else {
            BaseType::resize(NewSize);
        }
    }

    /**
     * @brief uBLAS-style clear
     * @details Value-initializes every entry; the size is kept.
     */
    void clear()
    {
        std::fill(begin(), end(), TDataType());
    }

    /**
     * @brief Swaps the contents with another vector
     * @details Delegates to Eigen's swap, which for dynamic vectors of the same type only exchanges the data
     * pointers.
     * @param rOther Vector whose contents are exchanged with this one.
     */
    void swap(EigenVector& rOther)
    {
        BaseType::swap(rOther);
    }

    /**
     * @brief Swaps the contents of two vectors
     * @details Non-member overload found by argument-dependent lookup; it forwards to the swap member.
     * @param rFirst Vector whose contents are exchanged with rSecond.
     * @param rSecond Vector whose contents are exchanged with rFirst.
     */
    friend void swap(EigenVector& rFirst, EigenVector& rSecond)
    {
        rFirst.swap(rSecond);
    }

    /**
     * @brief uBLAS-style element insertion
     * @details Assigns rValue to the entry at the given index (no entry is added) and returns a reference to it.
     * @param Index Position of the entry to assign.
     * @param rValue Value assigned to the entry.
     * @return Reference to the assigned entry.
     */
    reference insert_element(const std::size_t Index, const_reference rValue)
    {
        return ((*this)[Index] = rValue);
    }

    /**
     * @brief uBLAS-style element erasure
     * @details Value-initializes the entry at the given index (no entry is removed).
     * @param Index Position of the entry to reset.
     */
    void erase_element(const std::size_t Index)
    {
        (*this)[Index] = TDataType();
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief uBLAS-style size
     * @details Unsigned, so the pervasive comparisons and loops against std::size_t counters stay warning-free
     * (Eigen's own size() is signed).
     * @return Number of entries.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }

    /**
     * @brief uBLAS-style emptiness check
     * @details True when the vector has no entries.
     * @return True if the vector has no entries.
     */
    bool empty() const { return BaseType::size() == 0; }

    /**
     * @brief Contiguous storage view of the entries
     * @details Provides the uBLAS data() surface and an implicit conversion to the raw pointer.
     * @return Storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<TDataType> data() { return {BaseType::data(), size()}; }

    /**
     * @brief Contiguous read-only storage view of the entries
     * @details Const version of data(): the view refers to const entries.
     * @return Read-only storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<const TDataType> data() const { return {BaseType::data(), size()}; }

    /**
     * @brief Iterator to the first entry
     * @details Pointer to the contiguous storage.
     * @return Pointer to the first entry.
     */
    iterator begin() { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry
     * @details Const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator begin() const { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry, explicitly const
     * @details Same as the const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator cbegin() const { return BaseType::data(); }

    /**
     * @brief Iterator past the last entry
     * @details Pointer one past the end of the contiguous storage.
     * @return Pointer one past the last entry.
     */
    iterator end() { return BaseType::data() + BaseType::size(); }

    /**
     * @brief Read-only iterator past the last entry
     * @details Const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator end() const { return BaseType::data() + BaseType::size(); }

    /**
     * @brief Read-only iterator past the last entry, explicitly const
     * @details Same as the const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator cend() const { return BaseType::data() + BaseType::size(); }

    /**
     * @brief Reverse iterator to the last entry
     * @details Reverse iterator that starts at end().
     * @return Reverse iterator starting at the last entry.
     */
    reverse_iterator rbegin() { return reverse_iterator(end()); }

    /**
     * @brief Read-only reverse iterator to the last entry
     * @details Const version of rbegin().
     * @return Read-only reverse iterator starting at the last entry.
     */
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }

    /**
     * @brief Reverse iterator past the first entry
     * @details Reverse iterator that ends at begin().
     * @return Reverse iterator past the first entry.
     */
    reverse_iterator rend() { return reverse_iterator(begin()); }

    /**
     * @brief Read-only reverse iterator past the first entry
     * @details Const version of rend().
     * @return Read-only reverse iterator past the first entry.
     */
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    ///@}
};

/**
 * @brief Eigen-backed implementation of the Kratos fixed-size dense matrix (BoundedMatrix)
 * @details Mirrors the uBLAS bounded_matrix member surface (size1/size2, resize, operator(), the assignment protocol,
 * iterators) so generic Kratos code compiles unchanged. Unlike ublas::bounded_matrix the dimensions are strictly
 * fixed: resizing to anything but (TSize1, TSize2) is an error (use a dynamic Matrix where a variable size is
 * needed).
 * @tparam T Scalar type stored in the matrix.
 * @tparam TSize1 Fixed number of rows.
 * @tparam TSize2 Fixed number of columns.
 */
template<class T, std::size_t TSize1, std::size_t TSize2>
class EigenBoundedMatrix : public Eigen::Matrix<T, static_cast<int>(TSize1), static_cast<int>(TSize2), Internals::EigenBoundedMatrixOptions(TSize1, TSize2)>
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * @brief Pointer definition of EigenBoundedMatrix
     * @details Defines the Pointer, SharedPointer, WeakPointer and UniquePointer aliases of this class.
     */
    KRATOS_CLASS_POINTER_DEFINITION(EigenBoundedMatrix);

    /// Plain Eigen matrix type this class derives from. Fixed-size TSize1 x TSize2 Eigen matrix of T whose storage
    /// options (layout and no over-alignment) come from Internals::EigenBoundedMatrixOptions.
    using BaseType = Eigen::Matrix<T, static_cast<int>(TSize1), static_cast<int>(TSize2), Internals::EigenBoundedMatrixOptions(TSize1, TSize2)>;

    /// Unsigned type used for sizes and indices. Alias of std::size_t, as in the uBLAS bounded matrix.
    using size_type = std::size_t;

    /// Signed type used for differences between positions. Alias of std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;

    /// Scalar type of the stored entries. Alias of T, as in the uBLAS bounded matrix.
    using value_type = T;

    /// Read-only reference to a stored entry. Alias of const T&.
    using const_reference = const T&;

    /// Mutable reference to a stored entry. Alias of T&.
    using reference = T&;

    /// Pointer to a stored entry. Alias of T*.
    using pointer = T*;

    /// Iterator over the entries. Plain pointer, since the storage is contiguous (see begin() and end()).
    using iterator = T*;

    /// Read-only iterator over the entries. Plain pointer to const, since the storage is contiguous (see begin() and
    /// end()).
    using const_iterator = const T*;

    /// Type of the contiguous storage array exposed by data(). Internals::StorageView of the entries, a raw pointer
    /// plus the length that also converts implicitly to the raw pointer. It plays the role of the uBLAS storage array.
    using array_type = Internals::StorageView<T>;

    /**
     * @brief Static bound on the number of rows
     * @details Equal to TSize1. Part of the ublas::bounded_matrix interface, named by generic code (e.g. SVDUtils).
     */
    static constexpr size_type max_size1 = TSize1;

    /**
     * @brief Static bound on the number of columns
     * @details Equal to TSize2. Part of the ublas::bounded_matrix interface, named by generic code (e.g. SVDUtils).
     */
    static constexpr size_type max_size2 = TSize2;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates a matrix of the full (TSize1 x TSize2) size, values NOT initialized.
     */
    EigenBoundedMatrix() : BaseType() {}

    /**
     * @brief Constructor with the extents
     * @details The size arguments are only there for interface compatibility: they must equal (TSize1, TSize2), which
     * is checked when KRATOS_DEBUG is defined. The values are NOT initialized.
     * @param Size1 Requested number of rows; must equal TSize1.
     * @param Size2 Requested number of columns; must equal TSize2.
     */
    EigenBoundedMatrix(size_type Size1, size_type Size2) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size1 != TSize1 || Size2 != TSize2) << "Wrong size in the construction of a bounded matrix [ requested (" << Size1 << "," << Size2 << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
    }

    /**
     * @brief Copy constructor
     * @details Copies every entry of rOther.
     * @param rOther Matrix to copy.
     */
    EigenBoundedMatrix(const EigenBoundedMatrix& rOther) : BaseType(rOther) {}

    /**
     * @brief Constructor from any Eigen expression of compatible static size
     * @details The constraint keeps incompatible constructions out of overload sets. The cast mirrors the uBLAS
     * element-wise conversion semantics for expressions of a different scalar type.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression whose values, cast to T, initialize the matrix.
     */
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Entry access (const version)
     * @details Returns a const reference to the entry (I, J). The indices are checked against the static size when
     * KRATOS_DEBUG is defined.
     * @param I Row index.
     * @param J Column index.
     * @return Const reference to the entry at row I and column J.
     */
    const_reference operator()(size_type I, size_type J) const
    {
        KRATOS_DEBUG_ERROR_IF(I >= TSize1 || J >= TSize2) << "Index greater than the size of the bounded matrix - indices are (" << I << "," << J << ")" << std::endl;
        return BaseType::coeff(I, J);
    }

    /**
     * @brief Entry access
     * @details Returns a reference to the entry (I, J). The indices are checked against the static size when
     * KRATOS_DEBUG is defined.
     * @param I Row index.
     * @param J Column index.
     * @return Reference to the entry at row I and column J.
     */
    reference operator()(size_type I, size_type J)
    {
        KRATOS_DEBUG_ERROR_IF(I >= TSize1 || J >= TSize2) << "Index greater than the size of the bounded matrix - indices are (" << I << "," << J << ")" << std::endl;
        return BaseType::coeffRef(I, J);
    }

    /**
     * @brief Copy assignment operator
     * @details Copies every entry of rOther.
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    EigenBoundedMatrix& operator=(const EigenBoundedMatrix& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /**
     * @brief Assignment from any Eigen expression of compatible static size
     * @details The assignment evaluates through a temporary (.eval()), reproducing the alias-safe uBLAS operator=
     * semantics (e.g. M = trans(M)); the temporary-free fast path remains noalias(M) = expr, as in uBLAS.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to assign to this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>().eval());
        return *this;
    }

    /**
     * @brief Addition assignment from any Eigen expression of compatible static size
     * @details The expression is cast to T and evaluated through a temporary before being added, so aliasing with
     * this matrix is safe.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to add to this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator+=(rExpression.template cast<T>().eval());
        return *this;
    }

    /**
     * @brief Subtraction assignment from any Eigen expression of compatible static size
     * @details The expression is cast to T and evaluated through a temporary before being subtracted, so aliasing
     * with this matrix is safe.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to subtract from this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator-=(rExpression.template cast<T>().eval());
        return *this;
    }

    // Scalar operators keep the base class implementations (no aliasing there)
    using BaseType::operator*=;
    using BaseType::operator/=;

    /**
     * @brief Equality comparison
     * @details True when every entry equals the corresponding entry of rOther (exact comparison with std::equal).
     * @param rOther Matrix to compare with.
     * @return True if every entry equals the corresponding entry of rOther.
     */
    bool operator==(const EigenBoundedMatrix& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Interface-compatibility resize
     * @details The dimensions are fixed, so resizing to a different size is an error (unlike ublas::bounded_matrix,
     * which can carry runtime sizes below the static maximum). Resizing to (TSize1, TSize2) is accepted and keeps the
     * entries when Preserve is true; otherwise they are value-initialized.
     * @param Size1 Requested number of rows; must equal TSize1.
     * @param Size2 Requested number of columns; must equal TSize2.
     * @param Preserve If true (the default) the entries are kept, otherwise they are value-initialized.
     */
    void resize(size_type Size1, size_type Size2, bool Preserve = true)
    {
        KRATOS_ERROR_IF(Size1 != TSize1 || Size2 != TSize2) << "A bounded matrix of static size (" << TSize1 << "," << TSize2 << ") cannot be resized to (" << Size1 << "," << Size2 << ") under the Eigen backend." << std::endl;
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    /**
     * @brief uBLAS assignment protocol: assignment without a protective temporary
     * @details The expression is cast to T and assigned through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to assign to this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenBoundedMatrix& assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() = rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place addition without a protective temporary
     * @details The expression is cast to T and added to this matrix through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to add to this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenBoundedMatrix& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() += rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place subtraction without a protective temporary
     * @details The expression is cast to T and subtracted from this matrix through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to subtract from this matrix, cast to T.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenBoundedMatrix& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() -= rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief Swaps the contents with another matrix
     * @details Exchanges the entries one by one, since the storage is inline and no pointers can be swapped. Swapping
     * a matrix with itself does nothing.
     * @param rOther Matrix whose contents are exchanged with this one.
     */
    void swap(EigenBoundedMatrix& rOther)
    {
        if (this != &rOther) {
            std::swap_ranges(begin(), end(), rOther.begin());
        }
    }

    /**
     * @brief Swaps the contents of two matrices
     * @details Non-member overload found by argument-dependent lookup; it forwards to the swap member.
     * @param rFirst Matrix whose contents are exchanged with rSecond.
     * @param rSecond Matrix whose contents are exchanged with rFirst.
     */
    friend void swap(EigenBoundedMatrix& rFirst, EigenBoundedMatrix& rSecond)
    {
        rFirst.swap(rSecond);
    }

    /**
     * @brief uBLAS-style clear
     * @details Value-initializes every entry; the size is fixed.
     */
    void clear()
    {
        std::fill(begin(), end(), value_type());
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Number of rows
     * @details Always TSize1, known at compile time.
     * @return Number of rows, always TSize1.
     */
    constexpr size_type size1() const { return TSize1; }

    /**
     * @brief Number of columns
     * @details Always TSize2, known at compile time.
     * @return Number of columns, always TSize2.
     */
    constexpr size_type size2() const { return TSize2; }

    /**
     * @brief Contiguous storage view of the entries (row-major)
     * @details Provides the uBLAS data() surface and an implicit conversion to the raw pointer.
     * @return Storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<T> data() { return {BaseType::data(), TSize1 * TSize2}; }

    /**
     * @brief Contiguous read-only storage view of the entries (row-major)
     * @details Const version of data(): the view refers to const entries.
     * @return Read-only storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<const T> data() const { return {BaseType::data(), TSize1 * TSize2}; }

    /**
     * @brief Iterator to the first entry
     * @details Pointer to the contiguous storage.
     * @return Pointer to the first entry.
     */
    iterator begin() { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry
     * @details Const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator begin() const { return BaseType::data(); }

    /**
     * @brief Iterator past the last entry
     * @details Pointer one past the end of the contiguous storage (TSize1 * TSize2 entries after begin()).
     * @return Pointer one past the last entry.
     */
    iterator end() { return BaseType::data() + TSize1 * TSize2; }

    /**
     * @brief Read-only iterator past the last entry
     * @details Const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator end() const { return BaseType::data() + TSize1 * TSize2; }

    ///@}
}; // Class EigenBoundedMatrix

/**
 * @brief Eigen-backed implementation of the Kratos fixed-size dense vector alias (BoundedVector)
 * @details Same surface as the Eigen array_1d but a distinct type, since generic code (e.g. the serializer) overloads
 * on array_1d and BoundedVector simultaneously. Unlike ublas::bounded_vector the size is strictly fixed: resizing to
 * anything but TSize is an error.
 * @tparam T Scalar type stored in the vector.
 * @tparam TSize Fixed number of entries.
 */
template<class T, std::size_t TSize>
class EigenBoundedVector : public Eigen::Matrix<T, static_cast<int>(TSize), 1, Eigen::ColMajor | Eigen::DontAlign>
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * @brief Pointer definition of EigenBoundedVector
     * @details Defines the Pointer, SharedPointer, WeakPointer and UniquePointer aliases of this class.
     */
    KRATOS_CLASS_POINTER_DEFINITION(EigenBoundedVector);

    /// Plain Eigen vector type this class derives from. Fixed-size column vector of TSize entries of T, column-major
    /// and without over-alignment (DontAlign).
    using BaseType = Eigen::Matrix<T, static_cast<int>(TSize), 1, Eigen::ColMajor | Eigen::DontAlign>;

    /// Unsigned type used for sizes and indices. Alias of std::size_t, as in the uBLAS bounded vector.
    using size_type = std::size_t;

    /// Signed type used for differences between positions. Alias of std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;

    /// Scalar type of the stored entries. Alias of T, as in the uBLAS bounded vector.
    using value_type = T;

    /// Read-only reference to a stored entry. Alias of const T&.
    using const_reference = const T&;

    /// Mutable reference to a stored entry. Alias of T&.
    using reference = T&;

    /// Pointer to a stored entry. Alias of T*.
    using pointer = T*;

    /// Iterator over the entries. Plain pointer, since the storage is contiguous (see begin() and end()).
    using iterator = T*;

    /// Read-only iterator over the entries. Plain pointer to const, since the storage is contiguous (see begin() and
    /// end()).
    using const_iterator = const T*;

    /// Reverse iterator over the entries. std::reverse_iterator of iterator (see rbegin() and rend()).
    using reverse_iterator = std::reverse_iterator<iterator>;

    /// Read-only reverse iterator over the entries. std::reverse_iterator of const_iterator (see rbegin() and rend()).
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    /// Type of the contiguous storage array exposed by data(). Internals::StorageView of the entries, a raw pointer
    /// plus the length that also converts implicitly to the raw pointer. It plays the role of the uBLAS storage array.
    using array_type = Internals::StorageView<T>;

    /**
     * @brief Static bound on the number of entries
     * @details Equal to TSize. Part of the ublas::bounded_vector interface.
     */
    static constexpr size_type max_size = TSize;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates a vector of the full size, values NOT initialized.
     */
    EigenBoundedVector() : BaseType() {}

    /**
     * @brief Constructor with the size
     * @details The size argument is only there for interface compatibility: it must equal TSize, which is checked
     * when KRATOS_DEBUG is defined. The values are NOT initialized.
     * @param Size Requested number of entries; must equal TSize.
     */
    explicit EigenBoundedVector(size_type Size) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size != TSize) << "Wrong size in the construction of a bounded vector [ requested " << Size << ", bounded size " << TSize << " ]." << std::endl;
    }

    /**
     * @brief Fill constructor
     * @details Creates a vector of the full size and sets its first Size entries to rValue. Size must not exceed
     * TSize, which is checked when KRATOS_DEBUG is defined. The remaining entries are NOT initialized.
     * @param Size Number of leading entries set to rValue; must not exceed TSize.
     * @param rValue Value assigned to the first Size entries.
     */
    EigenBoundedVector(size_type Size, const value_type& rValue) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size > TSize) << "Given size is greater than the size of the bounded vector!" << std::endl;
        std::fill(begin(), begin() + Size, rValue);
    }

    /**
     * @brief Copy constructor
     * @details Copies every entry of rOther.
     * @param rOther Vector to copy.
     */
    EigenBoundedVector(const EigenBoundedVector& rOther) : BaseType(rOther) {}

    /**
     * @brief Constructor from any Eigen expression of compatible static size
     * @details The expression is cast to T.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression whose values, cast to T, initialize the vector.
     */
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Entry access (const version)
     * @details Returns a const reference to the entry i. The index is checked against TSize when KRATOS_DEBUG is
     * defined.
     * @param i Index of the entry.
     * @return Const reference to the entry i.
     */
    const_reference operator()(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    /**
     * @brief Entry access
     * @details Returns a reference to the entry i. The index is checked against TSize when KRATOS_DEBUG is defined.
     * @param i Index of the entry.
     * @return Reference to the entry i.
     */
    reference operator()(size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    /**
     * @brief Subscript access (const version)
     * @details Returns a const reference to the entry i. The index is checked against TSize when KRATOS_DEBUG is
     * defined.
     * @param i Index of the entry.
     * @return Const reference to the entry i.
     */
    const_reference operator[](size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    /**
     * @brief Subscript access
     * @details Returns a reference to the entry i. The index is checked against TSize when KRATOS_DEBUG is defined.
     * @param i Index of the entry.
     * @return Reference to the entry i.
     */
    reference operator[](size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    /**
     * @brief Copy assignment operator
     * @details Copies every entry of rOther.
     * @param rOther Vector to copy.
     * @return Reference to this vector.
     */
    EigenBoundedVector& operator=(const EigenBoundedVector& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /**
     * @brief Assignment from any Eigen expression of compatible static size
     * @details Alias-safe: the expression is cast to T and evaluated through a temporary, as the ublas operator=.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to assign to this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>().eval());
        return *this;
    }

    /**
     * @brief Addition assignment from any Eigen expression of compatible static size
     * @details The expression is cast to T and evaluated through a temporary before being added, so aliasing with
     * this vector is safe.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to add to this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator+=(rExpression.template cast<T>().eval());
        return *this;
    }

    /**
     * @brief Subtraction assignment from any Eigen expression of compatible static size
     * @details The expression is cast to T and evaluated through a temporary before being subtracted, so aliasing
     * with this vector is safe.
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to subtract from this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator-=(rExpression.template cast<T>().eval());
        return *this;
    }

    // Scalar operators keep the base class implementations (no aliasing there)
    using BaseType::operator*=;
    using BaseType::operator/=;

    /**
     * @brief Equality comparison
     * @details True when every entry equals the corresponding entry of rOther (exact comparison with std::equal).
     * @param rOther Vector to compare with.
     * @return True if every entry equals the corresponding entry of rOther.
     */
    bool operator==(const EigenBoundedVector& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Interface-compatibility resize
     * @details The size is fixed, so resizing to a different size is an error (unlike ublas::bounded_vector).
     * Resizing to TSize is accepted and keeps the entries when Preserve is true; otherwise they are
     * value-initialized.
     * @param Size Requested number of entries; must equal TSize.
     * @param Preserve If true (the default) the entries are kept, otherwise they are value-initialized.
     */
    void resize(size_type Size, bool Preserve = true)
    {
        KRATOS_ERROR_IF(Size != TSize) << "A bounded vector of static size " << TSize << " cannot be resized to " << Size << " under the Eigen backend." << std::endl;
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    /**
     * @brief uBLAS assignment protocol: assignment without a protective temporary
     * @details The expression is cast to T and assigned through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to assign to this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenBoundedVector& assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() = rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place addition without a protective temporary
     * @details The expression is cast to T and added to this vector through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to add to this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenBoundedVector& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() += rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief uBLAS assignment protocol: in-place subtraction without a protective temporary
     * @details The expression is cast to T and subtracted from this vector through noalias().
     * @tparam TDerived Eigen expression type.
     * @param rExpression Expression to subtract from this vector, cast to T.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenBoundedVector& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() -= rExpression.template cast<T>();
        return *this;
    }

    /**
     * @brief Swaps the contents with another vector
     * @details Exchanges the entries one by one, since the storage is inline and no pointers can be swapped. Swapping
     * a vector with itself does nothing.
     * @param rOther Vector whose contents are exchanged with this one.
     */
    void swap(EigenBoundedVector& rOther)
    {
        if (this != &rOther) {
            std::swap_ranges(begin(), end(), rOther.begin());
        }
    }

    /**
     * @brief Swaps the contents of two vectors
     * @details Non-member overload found by argument-dependent lookup; it forwards to the swap member.
     * @param rFirst Vector whose contents are exchanged with rSecond.
     * @param rSecond Vector whose contents are exchanged with rFirst.
     */
    friend void swap(EigenBoundedVector& rFirst, EigenBoundedVector& rSecond)
    {
        rFirst.swap(rSecond);
    }

    /**
     * @brief uBLAS-style clear
     * @details Value-initializes every entry; the size is fixed.
     */
    void clear()
    {
        std::fill(begin(), end(), value_type());
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Number of entries
     * @details Always TSize, known at compile time.
     * @return Number of entries, always TSize.
     */
    constexpr size_type size() const
    {
        return TSize;
    }

    /**
     * @brief Contiguous storage view of the entries
     * @details Provides the uBLAS data() surface and an implicit conversion to the raw pointer.
     * @return Storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<T> data() { return {BaseType::data(), TSize}; }

    /**
     * @brief Contiguous read-only storage view of the entries
     * @details Const version of data(): the view refers to const entries.
     * @return Read-only storage view of the entries: pointer to the first one and number of entries.
     */
    Internals::StorageView<const T> data() const { return {BaseType::data(), TSize}; }

    /**
     * @brief Iterator to the first entry
     * @details Pointer to the contiguous storage.
     * @return Pointer to the first entry.
     */
    iterator begin() { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry
     * @details Const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator begin() const { return BaseType::data(); }

    /**
     * @brief Read-only iterator to the first entry, explicitly const
     * @details Same as the const version of begin().
     * @return Read-only pointer to the first entry.
     */
    const_iterator cbegin() const { return BaseType::data(); }

    /**
     * @brief Iterator past the last entry
     * @details Pointer one past the end of the contiguous storage (TSize entries after begin()).
     * @return Pointer one past the last entry.
     */
    iterator end() { return BaseType::data() + TSize; }

    /**
     * @brief Read-only iterator past the last entry
     * @details Const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator end() const { return BaseType::data() + TSize; }

    /**
     * @brief Read-only iterator past the last entry, explicitly const
     * @details Same as the const version of end().
     * @return Read-only pointer one past the last entry.
     */
    const_iterator cend() const { return BaseType::data() + TSize; }

    /**
     * @brief Reverse iterator to the last entry
     * @details Reverse iterator that starts at end().
     * @return Reverse iterator starting at the last entry.
     */
    reverse_iterator rbegin() { return reverse_iterator(end()); }

    /**
     * @brief Read-only reverse iterator to the last entry
     * @details Const version of rbegin().
     * @return Read-only reverse iterator starting at the last entry.
     */
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }

    /**
     * @brief Reverse iterator past the first entry
     * @details Reverse iterator that ends at begin().
     * @return Reverse iterator past the first entry.
     */
    reverse_iterator rend() { return reverse_iterator(begin()); }

    /**
     * @brief Read-only reverse iterator past the first entry
     * @details Const version of rend().
     * @return Read-only reverse iterator past the first entry.
     */
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    ///@}
}; // Class EigenBoundedVector

///@}
///@name Lazy uBLAS-style factories
///@{
// zero_matrix, identity_matrix, scalar_matrix, zero_vector, scalar_vector and
// unit_vector are, as in uBLAS, lazy expressions: no storage is allocated and
// they can be assigned to (or combined with) any dense type of matching size:
//     Matrix A = ZeroMatrix(n, m);   noalias(A) = ZeroMatrix(n, m);
//     A += IdentityMatrix(n);        BoundedMatrix<double,3,3> B = ZeroMatrix(3, 3);
// They derive from Eigen's nullary expressions so the whole Eigen expression
// machinery applies to them; the classes only add the uBLAS constructors.

namespace Internals
{
/**
 * @brief Nullary functor of the unit vector e_i
 * @details Evaluates to one at the entry mIndex and to zero everywhere else. It is the functor of the Eigen nullary
 * expression behind unit_vector.
 * @tparam T Scalar type of the entries.
 */
template<class T>
struct UnitVectorOp
{
    /**
     * @brief Position of the unit entry
     * @details Index of the only entry of the vector that evaluates to one.
     */
    Eigen::Index mIndex;

    /**
     * @brief Value of the entry at the given index
     * @details Returns T(1) when Index equals mIndex and T(0) otherwise.
     * @param Index Position of the entry to evaluate.
     * @return One (T(1)) when Index equals mIndex, zero (T(0)) otherwise.
     */
    T operator()(const Eigen::Index Index) const { return Index == mIndex ? T(1) : T(0); }
};
} // namespace Internals

/**
 * @brief Lazy zero matrix (uBLAS zero_matrix)
 * @details zero_matrix(n) is n x n and zero_matrix(n, m) is n x m. No storage is allocated: it derives from an Eigen
 * nullary expression with the constant value zero, so it can be assigned to (or combined with) any dense type of
 * matching size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class zero_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Constant-value expression (scalar_constant_op, set to
    /// zero) with the shape of a dynamic row-major matrix.
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Matrix to copy.
     */
    zero_matrix(const zero_matrix& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    zero_matrix& operator=(const zero_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0 x 0) zero matrix.
     */
    zero_matrix() : BaseType(0, 0, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    /**
     * @brief Square constructor
     * @details Creates a Size x Size zero matrix.
     * @param Size Number of rows and of columns.
     */
    explicit zero_matrix(const std::size_t Size) : BaseType(Size, Size, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    /**
     * @brief Constructor with the extents
     * @details Creates a Size1 x Size2 zero matrix.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    zero_matrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    /**
     * @brief Number of rows
     * @details Returns the row count as an unsigned value.
     * @return Number of rows.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns
     * @details Returns the column count as an unsigned value.
     * @return Number of columns.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/**
 * @brief Lazy identity matrix (uBLAS identity_matrix)
 * @details identity_matrix(n) is n x n and identity_matrix(n, m) is n x m, with ones on the main diagonal and zeros
 * elsewhere. No storage is allocated: it derives from an Eigen nullary expression (scalar_identity_op), so it can be
 * assigned to (or combined with) any dense type of matching size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class identity_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Identity expression (scalar_identity_op) with the shape
    /// of a dynamic row-major matrix.
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<T>, Internals::DynamicRowMajorMatrix<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Matrix to copy.
     */
    identity_matrix(const identity_matrix& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    identity_matrix& operator=(const identity_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0 x 0) identity matrix.
     */
    identity_matrix() : BaseType(0, 0) {}

    /**
     * @brief Square constructor
     * @details Creates a Size x Size identity matrix.
     * @param Size Number of rows and of columns.
     */
    explicit identity_matrix(const std::size_t Size) : BaseType(Size, Size) {}

    /**
     * @brief Constructor with the extents
     * @details Creates a Size1 x Size2 identity matrix (ones on the main diagonal).
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    identity_matrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /**
     * @brief Number of rows
     * @details Returns the row count as an unsigned value.
     * @return Number of rows.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns
     * @details Returns the column count as an unsigned value.
     * @return Number of columns.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/**
 * @brief Lazy constant matrix (uBLAS scalar_matrix)
 * @details Every entry equals the given value (1 by default, as in uBLAS). No storage is allocated: it derives from an
 * Eigen nullary expression with a constant value, so it can be assigned to (or combined with) any dense type of
 * matching size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class scalar_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Constant-value expression (scalar_constant_op) with the
    /// shape of a dynamic row-major matrix.
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Matrix to copy.
     */
    scalar_matrix(const scalar_matrix& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    scalar_matrix& operator=(const scalar_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0 x 0) matrix whose value is 1.
     */
    scalar_matrix() : BaseType(0, 0, Eigen::internal::scalar_constant_op<T>(T(1))) {}

    /**
     * @brief Constructor with the extents and the value
     * @details Creates a Size1 x Size2 matrix whose entries all equal rValue (1 by default).
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     * @param rValue Value of every entry (1 by default).
     */
    scalar_matrix(const std::size_t Size1, const std::size_t Size2, const T& rValue = T(1)) : BaseType(Size1, Size2, Eigen::internal::scalar_constant_op<T>(rValue)) {}

    /**
     * @brief Number of rows
     * @details Returns the row count as an unsigned value.
     * @return Number of rows.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns
     * @details Returns the column count as an unsigned value.
     * @return Number of columns.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/**
 * @brief Lazy zero vector (uBLAS zero_vector)
 * @details zero_vector(n) is a column vector of n zeros. No storage is allocated: it derives from an Eigen nullary
 * expression with the constant value zero, so it can be assigned to (or combined with) any dense type of matching
 * size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class zero_vector : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Constant-value expression (scalar_constant_op, set to
    /// zero) with the shape of a dynamic column vector.
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Vector to copy.
     */
    zero_vector(const zero_vector& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Vector to copy.
     * @return Reference to this vector.
     */
    zero_vector& operator=(const zero_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0-length) zero vector.
     */
    zero_vector() : BaseType(0, 1, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    /**
     * @brief Constructor with the size
     * @details Creates a zero vector of Size entries.
     * @param Size Number of entries.
     */
    explicit zero_vector(const std::size_t Size) : BaseType(Size, 1, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    /**
     * @brief Number of entries
     * @details Returns the length as an unsigned value.
     * @return Number of entries.
     */
    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }
};

/**
 * @brief Lazy constant vector (uBLAS scalar_vector)
 * @details Every entry equals the given value (1 by default, as in uBLAS). No storage is allocated: it derives from an
 * Eigen nullary expression with a constant value, so it can be assigned to (or combined with) any dense type of
 * matching size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class scalar_vector : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Constant-value expression (scalar_constant_op) with the
    /// shape of a dynamic column vector.
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Vector to copy.
     */
    scalar_vector(const scalar_vector& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Vector to copy.
     * @return Reference to this vector.
     */
    scalar_vector& operator=(const scalar_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0-length) vector whose value is 1.
     */
    scalar_vector() : BaseType(0, 1, Eigen::internal::scalar_constant_op<T>(T(1))) {}

    /**
     * @brief Constructor with the size and the value
     * @details Creates a vector of Size entries that all equal rValue (1 by default).
     * @param Size Number of entries.
     * @param rValue Value of every entry (1 by default).
     */
    explicit scalar_vector(const std::size_t Size, const T& rValue = T(1)) : BaseType(Size, 1, Eigen::internal::scalar_constant_op<T>(rValue)) {}

    /**
     * @brief Number of entries
     * @details Returns the length as an unsigned value.
     * @return Number of entries.
     */
    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }
};

/**
 * @brief Lazy unit vector e_i (uBLAS unit_vector)
 * @details A column vector whose entry at the position given at construction is one and whose other entries are zero.
 * No storage is allocated: it derives from an Eigen nullary expression driven by Internals::UnitVectorOp, so it can
 * be assigned to (or combined with) any dense type of matching size.
 * @tparam T Scalar type of the entries (double by default).
 */
template<class T = double>
class unit_vector : public Eigen::CwiseNullaryOp<Internals::UnitVectorOp<T>, Internals::DynamicColumnVector<T>>
{
public:
    /// Eigen nullary expression type this class derives from. Expression driven by Internals::UnitVectorOp with the
    /// shape of a dynamic column vector.
    using BaseType = Eigen::CwiseNullaryOp<Internals::UnitVectorOp<T>, Internals::DynamicColumnVector<T>>;

    /// Scalar type of the entries. Alias of T.
    using value_type = T;

    /// Unsigned type used for sizes. Alias of std::size_t.
    using size_type = std::size_t;

    /**
     * @brief Copy constructor
     * @details Defaulted: copies the extents and the functor of rOther.
     * @param rOther Vector to copy.
     */
    unit_vector(const unit_vector& rOther) = default;

    /**
     * @brief Copy assignment operator
     * @details Assignable like the uBLAS factories: the Eigen nullary base has const extents and functor, so the value
     * object is destroyed and rebuilt in place from rOther (self-assignment does nothing).
     * @param rOther Vector to copy.
     * @return Reference to this vector.
     */
    unit_vector& operator=(const unit_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    /**
     * @brief Default constructor
     * @details Creates an empty (0-length) unit vector whose unit index is 0.
     */
    unit_vector() : BaseType(0, 1, Internals::UnitVectorOp<T>{0}) {}

    /**
     * @brief Constructor with the size and the position of the unit entry
     * @details Creates a vector of Size entries where the entry Index (0 by default) is one and the others are zero.
     * @param Size Number of entries.
     * @param Index Position of the entry equal to one (0 by default).
     */
    explicit unit_vector(const std::size_t Size, const std::size_t Index = 0) : BaseType(Size, 1, Internals::UnitVectorOp<T>{static_cast<Eigen::Index>(Index)}) {}

    /**
     * @brief Number of entries
     * @details Returns the length as an unsigned value.
     * @return Number of entries.
     */
    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Position of the unit entry
     * @details Returns the index stored in the nullary functor, that is the only entry equal to one.
     * @return Position of the entry equal to one.
     */
    std::size_t index() const { return static_cast<std::size_t>(this->functor().mIndex); }
};

///@}
///@name Input and output
///@{
// The Eigen-backed types print and parse in the boost::numeric::ublas text
// format ("[n](v0,v1,...)" and "[n1,n2]((a00,a01),(a10,a11))"): the format is
// part of the established Kratos IO surface (mdpa Matrix/Vector values, json
// checks, restart files, printed expectations in tests), so it must not
// change with the backend. These exact-type overloads win over Eigen's
// generic MatrixBase operator<<.

/**
 * @brief Writes an EigenVector in the uBLAS text format "[n](v0,v1,...)"
 * @details Prints the size in square brackets followed by the comma-separated entries in parentheses. The format is
 * part of the established Kratos IO surface (mdpa Matrix/Vector values, json checks, restart files, printed
 * expectations in tests), so it does not change with the backend. This exact-type overload wins over Eigen's generic
 * MatrixBase operator<<. Returns the output stream.
 * @tparam TDataType Scalar type stored in the vector.
 * @param rOStream Output stream to write to.
 * @param rV Vector to write.
 * @return Reference to the output stream rOStream.
 */
template<class TDataType>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenVector<TDataType>& rV)
{
    rOStream << '[' << rV.size() << "](";
    for (std::size_t i = 0; i < rV.size(); ++i) {
        if (i != 0) rOStream << ',';
        rOStream << rV[i];
    }
    rOStream << ')';
    return rOStream;
}

/**
 * @brief Reads an EigenVector from the uBLAS text format "[n](v0,v1,...)"
 * @details The failbit is set on any mismatch, as boost's own operator>> does. The vector is only assigned once the
 * whole input has been parsed, so it is left unchanged on failure. Returns the input stream.
 * @tparam TDataType Scalar type stored in the vector.
 * @param rIStream Input stream to read from.
 * @param rV Vector that receives the parsed entries; left unchanged on failure.
 * @return Reference to the input stream rIStream.
 */
template<class TDataType>
inline std::istream& operator>>(std::istream& rIStream, EigenVector<TDataType>& rV)
{
    char c;
    std::size_t size = 0;
    if (!(rIStream >> c) || c != '[') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    if (!(rIStream >> size)) return rIStream;
    if (!(rIStream >> c) || c != ']') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    if (!(rIStream >> c) || c != '(') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    EigenVector<TDataType> result(size);
    for (std::size_t i = 0; i < size; ++i) {
        if (i != 0 && (!(rIStream >> c) || c != ',')) { rIStream.setstate(std::ios_base::failbit); return rIStream; }
        if (!(rIStream >> result[i])) return rIStream;
    }
    if (!(rIStream >> c) || c != ')') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    rV = std::move(result);
    return rIStream;
}

/**
 * @brief Writes an EigenMatrix in the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))"
 * @details Prints the two sizes in square brackets followed by the comma-separated rows, each one a parenthesized
 * comma-separated list of entries. Like the vector one, this exact-type overload wins over Eigen's generic MatrixBase
 * operator<< and only covers the default row-major layout of EigenMatrix. Returns the output stream.
 * @tparam TDataType Scalar type stored in the matrix.
 * @param rOStream Output stream to write to.
 * @param rM Matrix to write.
 * @return Reference to the output stream rOStream.
 */
template<class TDataType>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenMatrix<TDataType>& rM)
{
    rOStream << '[' << rM.size1() << ',' << rM.size2() << "](";
    for (std::size_t i = 0; i < rM.size1(); ++i) {
        if (i != 0) rOStream << ',';
        rOStream << '(';
        for (std::size_t j = 0; j < rM.size2(); ++j) {
            if (j != 0) rOStream << ',';
            rOStream << rM(i, j);
        }
        rOStream << ')';
    }
    rOStream << ')';
    return rOStream;
}

/**
 * @brief Reads an EigenMatrix from the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))"
 * @details The failbit is set on any mismatch, as boost's own operator>> does. The matrix is only assigned once the
 * whole input has been parsed, so it is left unchanged on failure. Only the default row-major layout of EigenMatrix
 * is covered. Returns the input stream.
 * @tparam TDataType Scalar type stored in the matrix.
 * @param rIStream Input stream to read from.
 * @param rM Matrix that receives the parsed entries; left unchanged on failure.
 * @return Reference to the input stream rIStream.
 */
template<class TDataType>
inline std::istream& operator>>(std::istream& rIStream, EigenMatrix<TDataType>& rM)
{
    char c;
    std::size_t size1 = 0, size2 = 0;
    if (!(rIStream >> c) || c != '[') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    if (!(rIStream >> size1)) return rIStream;
    if (!(rIStream >> c) || c != ',') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    if (!(rIStream >> size2)) return rIStream;
    if (!(rIStream >> c) || c != ']') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    if (!(rIStream >> c) || c != '(') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    EigenMatrix<TDataType> result(size1, size2);
    for (std::size_t i = 0; i < size1; ++i) {
        if (i != 0 && (!(rIStream >> c) || c != ',')) { rIStream.setstate(std::ios_base::failbit); return rIStream; }
        if (!(rIStream >> c) || c != '(') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
        for (std::size_t j = 0; j < size2; ++j) {
            if (j != 0 && (!(rIStream >> c) || c != ',')) { rIStream.setstate(std::ios_base::failbit); return rIStream; }
            if (!(rIStream >> result(i, j))) return rIStream;
        }
        if (!(rIStream >> c) || c != ')') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    }
    if (!(rIStream >> c) || c != ')') { rIStream.setstate(std::ios_base::failbit); return rIStream; }
    rM = std::move(result);
    return rIStream;
}

/**
 * @brief Writes an EigenBoundedMatrix in the uBLAS text format "[N1,N2]((...),(...))"
 * @details Prints the static sizes in square brackets followed by the comma-separated rows, each one a parenthesized
 * comma-separated list of entries. Returns the output stream.
 * @tparam T Scalar type stored in the matrix.
 * @tparam TSize1 Fixed number of rows.
 * @tparam TSize2 Fixed number of columns.
 * @param rOStream Output stream to write to.
 * @param rMatrix Bounded matrix to write.
 * @return Reference to the output stream rOStream.
 */
template<class T, std::size_t TSize1, std::size_t TSize2>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenBoundedMatrix<T, TSize1, TSize2>& rMatrix)
{
    rOStream << '[' << TSize1 << ',' << TSize2 << "](";
    for (std::size_t i = 0; i < TSize1; ++i) {
        if (i > 0) rOStream << ',';
        rOStream << '(';
        for (std::size_t j = 0; j < TSize2; ++j) {
            if (j > 0) rOStream << ',';
            rOStream << rMatrix(i, j);
        }
        rOStream << ')';
    }
    rOStream << ')';
    return rOStream;
}

/**
 * @brief Writes an EigenBoundedVector in the uBLAS text format "[N](v0,v1,...)"
 * @details Prints the static size in square brackets followed by the comma-separated entries in parentheses. Returns
 * the output stream.
 * @tparam T Scalar type stored in the vector.
 * @tparam TSize Fixed number of entries.
 * @param rOStream Output stream to write to.
 * @param rVector Bounded vector to write.
 * @return Reference to the output stream rOStream.
 */
template<class T, std::size_t TSize>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenBoundedVector<T, TSize>& rVector)
{
    rOStream << '[' << TSize << "](";
    for (std::size_t i = 0; i < TSize; ++i) {
        if (i > 0) rOStream << ',';
        rOStream << rVector[i];
    }
    rOStream << ')';
    return rOStream;
}

///@}

// Contiguity contract (see the array_1d layout notes): the bounded types are
// plain T[N] blocks without over-alignment, so raw-buffer interoperability
// with the uBLAS storage is preserved.
static_assert(sizeof(EigenBoundedMatrix<double, 3, 3>) == 9 * sizeof(double), "EigenBoundedMatrix must be layout-compatible with double[Size1*Size2]");
static_assert(sizeof(EigenBoundedVector<double, 4>) == 4 * sizeof(double), "EigenBoundedVector must be layout-compatible with double[Size]");
static_assert(alignof(EigenBoundedMatrix<double, 2, 2>) == alignof(double), "The bounded types must not be over-aligned");

} // namespace Kratos
