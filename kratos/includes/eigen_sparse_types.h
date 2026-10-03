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
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "includes/eigen_dense_types.h"

namespace Kratos
{

/**
 * @brief Index type used for the Eigen sparse matrices.
 * @details Eigen requires a signed StorageIndex. The default is a 32-bit int:
 * sparse kernels are memory-bound (SpMV traffic drops from 16 to 12 bytes per
 * nonzero versus 64-bit indices), and 2^31 - 1 nonzeros per matrix (~16 GB of
 * values alone) comfortably covers shared-memory problems. For larger systems
 * configure with KRATOS_EIGEN_64BIT_INDICES=ON (or override the type directly
 * with -DKRATOS_EIGEN_INDEX_TYPE=<type>).
 */
#if defined(KRATOS_EIGEN_INDEX_TYPE)
using KratosEigenIndexType = KRATOS_EIGEN_INDEX_TYPE;
#elif defined(KRATOS_EIGEN_64BIT_INDICES)
using KratosEigenIndexType = std::ptrdiff_t;
#else
using KratosEigenIndexType = int;
#endif
static_assert(std::is_signed_v<KratosEigenIndexType>, "Eigen requires a signed sparse StorageIndex.");

///@name Kratos Classes
///@{

/**
 * @class EigenCompressedMatrix
 * @brief Row-major (CSR) Eigen sparse matrix with the uBLAS compressed_matrix member surface.
 * @details The compressed storage of Eigen::SparseMatrix in row-major mode is
 * exactly the CSR triplet of arrays that boost::numeric::ublas::compressed_matrix
 * exposes through index1_data() (row pointers), index2_data() (column indices)
 * and value_data() (values). This wrapper adds those accessors plus the
 * (rows, cols, nnz) constructor and set_filled() so the graph-construction and
 * assembly code of the builder-and-solvers works unchanged on Eigen storage.
 * @tparam TDataType The scalar type stored in the matrix (e.g. double).
 * @tparam TIndexType The signed index type used for the CSR storage (defaults to KratosEigenIndexType).
 * @warning As for the uBLAS type, element insertion through operator() on a
 * missing entry is an O(nnz-in-row) slow path and must not be used for assembly.
 */
template<class TDataType, class TIndexType = KratosEigenIndexType>
class EigenCompressedMatrix : public Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>;
    using value_type = TDataType;
    using size_type = std::size_t;
    // uBLAS-style storage-array typedefs (generic code uses e.g.
    // typename TMatrix::index_array_type::value_type to name the index type)
    using index_array_type = Internals::StorageView<TIndexType>;
    using value_array_type = Internals::StorageView<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: an empty (0 x 0) matrix.
    EigenCompressedMatrix() = default;

    /// Copy constructor: the copy gets the complete row pointers even when
    /// the source is still being built by ordered insertion.
    EigenCompressedMatrix(const EigenCompressedMatrix& rOther) : BaseType(rOther), mAppendRow(rOther.mAppendRow)
    {
        CompleteAppend();
    }

    EigenCompressedMatrix(EigenCompressedMatrix&& rOther) noexcept : BaseType(std::move(static_cast<BaseType&>(rOther))), mAppendRow(rOther.mAppendRow)
    {
        rOther.mAppendRow = -1;
    }

    /// Copy assignment (the target gets the complete row pointers, see the copy constructor).
    EigenCompressedMatrix& operator=(const EigenCompressedMatrix& rOther)
    {
        BaseType::operator=(rOther);
        mAppendRow = rOther.mAppendRow;
        CompleteAppend();
        return *this;
    }

    EigenCompressedMatrix& operator=(EigenCompressedMatrix&& rOther) noexcept
    {
        // Eigen moves by swapping the storage, so the build state goes along
        BaseType::operator=(std::move(static_cast<BaseType&>(rOther)));
        std::swap(mAppendRow, rOther.mAppendRow);
        return *this;
    }

    /// Swaps the storage and the build state with another matrix.
    void swap(EigenCompressedMatrix& rOther)
    {
        BaseType::swap(rOther);
        std::swap(mAppendRow, rOther.mAppendRow);
    }

    /// Allocates an empty Size1 x Size2 sparse matrix (no nonzeros reserved).
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /// uBLAS-style (rows, cols, nnz) constructor: allocates the compressed
    /// storage up front so the CSR arrays can be written directly.
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2, const std::size_t NNZ)
        : BaseType(Size1, Size2)
    {
        this->resizeNonZeros(NNZ);
    }

    /// Construction from any Eigen sparse expression.
    template<class TDerived>
    EigenCompressedMatrix(const Eigen::SparseMatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /// Construction from a dense Eigen expression: the nonzero entries are
    /// gathered into compressed storage.
    template<class TDerived>
    explicit EigenCompressedMatrix(const Eigen::MatrixBase<TDerived>& rExpression)
        : BaseType(rExpression.derived().sparseView())
    {
        this->makeCompressed();
    }

    /// Construction from a lazy zero matrix: the uBLAS idiom
    /// `SparseMatrixType m = ZeroMatrix(n, n);` only sets the dimensions.
    EigenCompressedMatrix(const zero_matrix<TDataType>& rZero) : BaseType(rZero.size1(), rZero.size2()) {}

    /// Construction from a lazy identity matrix: the uBLAS idiom
    /// `SparseMatrixType m = IdentityMatrix(n, n);` stores the unit diagonal.
    EigenCompressedMatrix(const identity_matrix<TDataType>& rIdentity) : BaseType(rIdentity.size1(), rIdentity.size2())
    {
        const std::size_t diagonal_size = std::min(size1(), size2());
        this->reserve(diagonal_size);
        for (std::size_t i = 0; i < diagonal_size; ++i) {
            push_back(i, i, TDataType(1));
        }
        CompleteAppend();
    }

    /// Element-wise conversion from a compressed matrix of another scalar or index type.
    template<class TOtherDataType, class TOtherIndexType>
    explicit EigenCompressedMatrix(const EigenCompressedMatrix<TOtherDataType, TOtherIndexType>& rOther)
        : BaseType(rOther.template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment from any Eigen sparse expression.
    template<class TDerived>
    EigenCompressedMatrix& operator=(const Eigen::SparseMatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        mAppendRow = -1;
        return *this;
    }

    /**
     * @brief Reference to a stored entry (ublas sparse_matrix_element).
     * @details Returned by the non-const operator(), so the uBLAS idioms
     * `A(i, j) = v`, `A(i, j) += v`, `double& r = A(i, j)` and
     * `AtomicAdd(A(i, j).ref(), v)` compile unchanged.
     */
    class ElementReference
    {
    public:
        explicit ElementReference(TDataType& rEntry) : mrEntry(rEntry) {}

        /// The referenced storage entry (as ublas sparse_matrix_element::ref()).
        TDataType& ref() const { return mrEntry; }

        operator TDataType&() const { return mrEntry; }

        ElementReference& operator=(const TDataType& rValue) { mrEntry = rValue; return *this; }
        ElementReference& operator=(const ElementReference& rOther) { mrEntry = rOther.mrEntry; return *this; }
        ElementReference& operator+=(const TDataType& rValue) { mrEntry += rValue; return *this; }
        ElementReference& operator-=(const TDataType& rValue) { mrEntry -= rValue; return *this; }
        ElementReference& operator*=(const TDataType& rValue) { mrEntry *= rValue; return *this; }
        ElementReference& operator/=(const TDataType& rValue) { mrEntry /= rValue; return *this; }

    private:
        TDataType& mrEntry;
    };

    /// uBLAS-style element access. As for ublas::compressed_matrix, a new
    /// entry of the last filled row (or of a later one) only shifts the
    /// entries of that row, so a row-by-row construction is O(1) (amortized)
    /// per entry when the columns come in increasing order; any other new
    /// entry is the O(nnz) slow path. The storage stays a packed CSR in both
    /// cases (Eigen's own insertion would leave it in uncompressed mode).
    ElementReference operator()(const std::size_t I, const std::size_t J)
    {
        if (IsTailRow(I, J)) {
            return ElementReference(TailRowEntry(I, J));
        }
        CompleteAppend();
        TDataType& r_entry = this->coeffRef(I, J);
        if (!this->isCompressed()) {
            this->makeCompressed();
            return ElementReference(this->coeffRef(I, J));
        }
        return ElementReference(r_entry);
    }

    /// uBLAS-style element access (const version).
    TDataType operator()(const std::size_t I, const std::size_t J) const { return this->coeff(I, J); }

    ///@}
    ///@name Access
    ///@{

    /// Number of rows.
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /// Number of columns.
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /// Number of stored entries (the *filled* count, as in ublas::compressed_matrix).
    std::size_t nnz() const { return static_cast<std::size_t>(this->nonZeros()); }

    /// Allocated storage capacity for nonzero entries (ublas nnz_capacity()).
    std::size_t nnz_capacity() const { return static_cast<std::size_t>(this->data().size()); }

    // The storage views span the *allocated* capacity (as the uBLAS storage
    // arrays do), so they can be written before the row pointers declare the
    // filled size.

    /// View over the CSR value array.
    auto value_data() { return Internals::StorageView<TDataType>(this->valuePtr(), this->data().size()); }
    /// View over the CSR row-pointer array (size1() + 1 entries).
    auto index1_data() { return Internals::StorageView<TIndexType>(this->outerIndexPtr(), size1() + 1); }
    /// View over the CSR column-index array.
    auto index2_data() { return Internals::StorageView<TIndexType>(this->innerIndexPtr(), this->data().size()); }
    /// View over the CSR value array (const version).
    auto value_data() const { return Internals::StorageView<const TDataType>(this->valuePtr(), this->data().size()); }
    /// View over the CSR row-pointer array (const version).
    auto index1_data() const { return Internals::StorageView<const TIndexType>(this->outerIndexPtr(), size1() + 1); }
    /// View over the CSR column-index array (const version).
    auto index2_data() const { return Internals::StorageView<const TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /// uBLAS-style filled row-pointer count (size1() + 1 once assembled).
    std::size_t filled1() const { return size1() + 1; }

    /// uBLAS-style filled nonzero count (nnz()).
    std::size_t filled2() const { return nnz(); }

    ///@}
    ///@name Iterators
    ///@{
    // Row-major traversal with the boost::numeric::ublas::compressed_matrix
    // iterator concept, so the ublas idiom compiles unchanged:
    //     for (auto i1 = m.begin1(); i1 != m.end1(); ++i1)
    //         for (auto i2 = i1.begin(); i2 != i1.end(); ++i2)
    //             ... i2.index1(), i2.index2(), *i2 ...
    // Only the row-major direction is provided (begin2()/iterator2 traversal
    // over a column has no users and is not meaningful for a CSR matrix).

    /**
     * @brief Iterator over the stored entries of a single row (ublas iterator2).
     * @details Walks the CSR arrays of one row; index1() is the row, index2()
     * the column of the entry it currently points at.
     * @tparam TIsConst Whether the iterator grants write access to the values.
     */
    template<bool TIsConst>
    class EntryIterator
    {
    public:
        using MatrixPointerType = std::conditional_t<TIsConst, const EigenCompressedMatrix*, EigenCompressedMatrix*>;
        using reference = std::conditional_t<TIsConst, const TDataType&, TDataType&>;
        using value_type = TDataType;

        EntryIterator(MatrixPointerType pMatrix, const std::size_t Row, const std::size_t Position)
            : mpMatrix(pMatrix), mRow(Row), mPosition(Position) {}

        /// Row of the entry currently pointed at.
        std::size_t index1() const { return mRow; }

        /// Column of the entry currently pointed at.
        std::size_t index2() const { return static_cast<std::size_t>(mpMatrix->innerIndexPtr()[mPosition]); }

        reference operator*() const { return mpMatrix->valuePtr()[mPosition]; }

        EntryIterator& operator++() { ++mPosition; return *this; }
        EntryIterator operator++(int) { EntryIterator copy(*this); ++mPosition; return copy; }

        bool operator==(const EntryIterator& rOther) const { return mPosition == rOther.mPosition; }
        bool operator!=(const EntryIterator& rOther) const { return !(*this == rOther); }

    private:
        MatrixPointerType mpMatrix;
        std::size_t mRow;
        std::size_t mPosition; /// Index into the CSR value/column arrays.
    };

    /**
     * @brief Iterator over the rows of the matrix (ublas iterator1).
     * @details begin()/end() yield the entry iterators of the current row.
     * @tparam TIsConst Whether the entry iterators grant write access.
     */
    template<bool TIsConst>
    class RowIterator
    {
    public:
        using MatrixPointerType = std::conditional_t<TIsConst, const EigenCompressedMatrix*, EigenCompressedMatrix*>;
        using EntryIteratorType = EntryIterator<TIsConst>;

        RowIterator(MatrixPointerType pMatrix, const std::size_t Row)
            : mpMatrix(pMatrix), mRow(Row) {}

        /// Index of the row currently pointed at.
        std::size_t index1() const { return mRow; }

        /// First stored entry of this row.
        EntryIteratorType begin() const
        {
            return EntryIteratorType(mpMatrix, mRow, static_cast<std::size_t>(mpMatrix->outerIndexPtr()[mRow]));
        }

        /// Past-the-last stored entry of this row.
        EntryIteratorType end() const
        {
            return EntryIteratorType(mpMatrix, mRow, static_cast<std::size_t>(mpMatrix->outerIndexPtr()[mRow + 1]));
        }

        RowIterator& operator++() { ++mRow; return *this; }
        RowIterator operator++(int) { RowIterator copy(*this); ++mRow; return copy; }

        bool operator==(const RowIterator& rOther) const { return mRow == rOther.mRow; }
        bool operator!=(const RowIterator& rOther) const { return !(*this == rOther); }

    private:
        MatrixPointerType mpMatrix;
        std::size_t mRow;
    };

    using iterator1 = RowIterator<false>;
    using const_iterator1 = RowIterator<true>;
    using iterator2 = EntryIterator<false>;
    using const_iterator2 = EntryIterator<true>;

    /// First row (the storage is packed first, so the CSR arrays the iterators
    /// walk are the compressed ones; a no-op when it already is compressed).
    iterator1 begin1()
    {
        complete_index1_data();
        return iterator1(this, 0);
    }

    /// Past-the-last row.
    iterator1 end1() { return iterator1(this, size1()); }

    /// First row (const version).
    const_iterator1 begin1() const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(this->isCompressed())
            << "Iterating a compressed matrix whose storage is not packed; call complete_index1_data() first." << std::endl;
        return const_iterator1(this, 0);
    }

    /// Past-the-last row (const version).
    const_iterator1 end1() const { return const_iterator1(this, size1()); }

    ///@}
    ///@name Operations
    ///@{

    /// uBLAS-style finalization after writing the CSR arrays by hand: the
    /// storage was already sized and the filled count is implied by the
    /// written row pointers, so this only validates consistency.
    void set_filled(const std::size_t FilledSize1, const std::size_t FilledNNZ)
    {
        CompleteAppend();
        KRATOS_DEBUG_ERROR_IF(FilledSize1 != size1() + 1 || FilledNNZ != nnz() || FilledNNZ > static_cast<std::size_t>(this->data().size()))
            << "set_filled(" << FilledSize1 << ", " << FilledNNZ
            << ") is inconsistent with the written compressed storage ("
            << size1() + 1 << ", " << nnz() << " of " << this->data().size()
            << " allocated)." << std::endl;
    }

    /// uBLAS-style resize; as in ublas::compressed_matrix the default
    /// preserves, while resize(m, n, false) discards values and structure.
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        CompleteAppend();
        if (Preserve) {
            this->conservativeResize(NewSize1, NewSize2);
        } else {
            BaseType::resize(NewSize1, NewSize2);
        }
    }

    // Brings in BaseType's reserve() overloads (single count, or per-inner-vector
    // counts) alongside the uBLAS-style two-argument one added below.
    using BaseType::reserve;

    /// uBLAS-style reservation of nonzero storage capacity: Preserve=false
    /// discards existing values and structure while allocating room for NNZ
    /// nonzeros; Preserve=true keeps the existing structure.
    void reserve(const std::size_t NNZ, const bool Preserve = true)
    {
        CompleteAppend();
        if (Preserve) {
            BaseType::reserve(static_cast<typename BaseType::Index>(NNZ));
        } else {
            this->resizeNonZeros(static_cast<typename BaseType::Index>(NNZ));
        }
    }

    /// uBLAS-style clear: removes all stored entries (the dimensions are kept).
    void clear()
    {
        mAppendRow = -1;
        this->setZero();
    }

    /// uBLAS-style ordered insertion, as ublas::compressed_matrix::push_back():
    /// entries given row by row with increasing columns are appended in O(1)
    /// (amortized).
    void push_back(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        (*this)(I, J) = Value;
    }

    /// uBLAS-style element insertion, as ublas::compressed_matrix::insert_element.
    void insert_element(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        (*this)(I, J) = Value;
    }

    /// uBLAS-style finalization, as ublas::compressed_matrix::complete_index1_data():
    /// sets the row pointers of the trailing rows left empty by the ordered
    /// insertion, and packs the storage after Eigen-native insertions
    /// (insert(), reserve()) left it uncompressed.
    void complete_index1_data()
    {
        CompleteAppend();
        this->makeCompressed();
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Row of the last entry appended by the ordered insertion while the row
    /// pointers of the rows after it are not yet written (-1 otherwise). As
    /// for ublas::compressed_matrix::filled1_, the row pointers past it are
    /// completed by the next non-ordered access or complete_index1_data().
    std::ptrdiff_t mAppendRow = -1;

    ///@}
    ///@name Private Operations
    ///@{

    /// Number of entries written so far (the packed filled count).
    std::size_t FilledEntries() const
    {
        const TIndexType* p_outer = this->outerIndexPtr();
        return static_cast<std::size_t>(mAppendRow >= 0 ? p_outer[mAppendRow + 1] : p_outer[size1()]);
    }

    /// Whether row I is the last row with stored entries or comes after it,
    /// so an entry of it can be inserted by shifting that row only.
    bool IsTailRow(const std::size_t I, const std::size_t J) const
    {
        if (!this->isCompressed() || I >= size1() || J >= size2()) {
            return false;
        }
        if (mAppendRow >= 0) {
            return I >= static_cast<std::size_t>(mAppendRow);
        }
        // All rows after I must be empty
        return static_cast<std::size_t>(this->outerIndexPtr()[I + 1]) == FilledEntries();
    }

    /// Entry (I, J) of the tail row I (IsTailRow(I, J) must hold), inserted
    /// with a zero value if not stored yet. The columns of the row are kept
    /// sorted, shifting the entries of that row only: an ordered insertion is
    /// an O(1) (amortized) append, as for ublas::compressed_matrix.
    TDataType& TailRowEntry(const std::size_t I, const std::size_t J)
    {
        const std::size_t filled = FilledEntries();
        TIndexType* p_outer = this->outerIndexPtr();
        const bool is_new_row = mAppendRow >= 0 && I > static_cast<std::size_t>(mAppendRow);
        const std::size_t row_begin = is_new_row ? filled : static_cast<std::size_t>(p_outer[I]);

        // Position of column J among the (sorted) columns of the row
        const TIndexType column = static_cast<TIndexType>(J);
        const TIndexType* p_row_begin = this->innerIndexPtr() + row_begin;
        const TIndexType* p_row_end = this->innerIndexPtr() + filled;
        const std::size_t position = row_begin + static_cast<std::size_t>(std::lower_bound(p_row_begin, p_row_end, column) - p_row_begin);
        if (position < filled && this->innerIndexPtr()[position] == column) {
            return this->data().value(position);
        }

        // The rows skipped since the last appended one are empty (when the
        // row pointers are complete, the rows up to I are already correct)
        if (mAppendRow >= 0) {
            for (std::size_t k = static_cast<std::size_t>(mAppendRow + 2); k <= I; ++k) {
                p_outer[k] = static_cast<TIndexType>(filled);
            }
        }
        // Reuse storage allocated up front ((rows, cols, nnz) constructor or
        // reserve(nnz, false)), growing it geometrically otherwise
        if (filled == static_cast<std::size_t>(this->data().size())) {
            this->data().append(TDataType(), column);
            p_outer = this->outerIndexPtr();
        }
        TDataType* p_values = this->valuePtr();
        TIndexType* p_columns = this->innerIndexPtr();
        std::copy_backward(p_values + position, p_values + filled, p_values + filled + 1);
        std::copy_backward(p_columns + position, p_columns + filled, p_columns + filled + 1);
        p_values[position] = TDataType();
        p_columns[position] = column;
        p_outer[I + 1] = static_cast<TIndexType>(filled + 1);
        p_outer[size1()] = static_cast<TIndexType>(filled + 1);
        // Once the last row is reached no row pointer is pending, so the
        // matrix is left complete (later concurrent reads of existing entries
        // then never write the build state)
        mAppendRow = (I + 1 == size1()) ? -1 : static_cast<std::ptrdiff_t>(I);
        return p_values[position];
    }

    /// Writes the row pointers of the rows after the last appended one.
    void CompleteAppend()
    {
        if (mAppendRow >= 0) {
            TIndexType* p_outer = this->outerIndexPtr();
            const TIndexType filled = p_outer[mAppendRow + 1];
            for (std::size_t k = static_cast<std::size_t>(mAppendRow + 2); k <= size1(); ++k) {
                p_outer[k] = filled;
            }
            mAppendRow = -1;
        }
    }

    ///@}
};

///@}

} // namespace Kratos
