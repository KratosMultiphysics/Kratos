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
 * @brief Index type used for the Eigen sparse matrices
 * @details Eigen requires a signed StorageIndex. The default (CMake option
 * KRATOS_EIGEN_64BIT_INDICES=ON) is a 64-bit std::ptrdiff_t, which does not limit
 * the number of nonzeros per matrix. A 32-bit int (KRATOS_EIGEN_64BIT_INDICES=OFF)
 * lowers the memory traffic of the memory-bound sparse kernels (SpMV traffic drops
 * from 16 to 12 bytes per nonzero), but caps the matrix at 2^31 - 1 nonzeros
 * (~16 GB of values alone). The type can also be overridden directly with
 * `-DKRATOS_EIGEN_INDEX_TYPE=<type>`.
 */
#if defined(KRATOS_EIGEN_INDEX_TYPE)
/// Index type taken from the KRATOS_EIGEN_INDEX_TYPE macro. Active when KRATOS_EIGEN_INDEX_TYPE is defined, which takes
/// precedence over KRATOS_EIGEN_64BIT_INDICES. The type must be signed, as checked by the static_assert below.
using KratosEigenIndexType = KRATOS_EIGEN_INDEX_TYPE;
#elif defined(KRATOS_EIGEN_64BIT_INDICES)
/// 64-bit index type. Active when KRATOS_EIGEN_64BIT_INDICES is defined (CMake option KRATOS_EIGEN_64BIT_INDICES=ON)
/// and KRATOS_EIGEN_INDEX_TYPE is not. std::ptrdiff_t does not limit the number of nonzeros per matrix.
using KratosEigenIndexType = std::ptrdiff_t;
#else
/// 32-bit index type. Active when neither KRATOS_EIGEN_INDEX_TYPE nor KRATOS_EIGEN_64BIT_INDICES is defined
/// (KRATOS_EIGEN_64BIT_INDICES=OFF). int lowers the memory traffic of the sparse kernels but caps the matrix at 2^31 -
/// 1 nonzeros.
using KratosEigenIndexType = int;
#endif
static_assert(std::is_signed_v<KratosEigenIndexType>, "Eigen requires a signed sparse StorageIndex.");

///@name Kratos Classes
///@{

/**
 * @brief Row-major (CSR) Eigen sparse matrix with the uBLAS compressed_matrix member surface
 * @details The compressed storage of Eigen::SparseMatrix in row-major mode is
 * exactly the CSR triplet of arrays that boost::numeric::ublas::compressed_matrix
 * exposes through index1_data() (row pointers), index2_data() (column indices)
 * and value_data() (values). This wrapper adds those accessors plus the
 * (rows, cols, nnz) constructor and set_filled() so the graph-construction and
 * assembly code of the builder-and-solvers works unchanged on Eigen storage.
 * Warning: as for the uBLAS type, element insertion through operator() on a missing
 * entry is an O(nnz) slow path and must not be used for assembly.
 * @tparam TDataType Scalar type stored in the matrix (e.g. double).
 * @tparam TIndexType Signed index type used for the CSR storage (defaults to KratosEigenIndexType).
 */
template<class TDataType, class TIndexType = KratosEigenIndexType>
class EigenCompressedMatrix : public Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Eigen sparse matrix this class derives from. Row-major (CSR) Eigen::SparseMatrix of TDataType scalars with
    /// TIndexType as storage index.
    using BaseType = Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>;

    /// Scalar type stored in the matrix. Alias of TDataType, as the uBLAS containers expose it.
    using value_type = TDataType;

    /// Type used for sizes and indices. Alias of std::size_t, the type returned by size1(), size2() and nnz() and taken
    /// by operator().
    using size_type = std::size_t;

    // uBLAS-style storage-array typedefs (generic code uses e.g.
    // typename TMatrix::index_array_type::value_type to name the index type)
    /// Type of the CSR index arrays. Kratos::span over TIndexType, the type returned by the non-const
    /// index1_data() and index2_data(). Generic code uses index_array_type::value_type to name the index type of the
    /// matrix.
    using index_array_type = Kratos::span<TIndexType>;

    /// Type of the CSR value array. Kratos::span over TDataType, the type returned by the non-const
    /// value_data().
    using value_array_type = Kratos::span<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Creates an empty (0 x 0) matrix.
     */
    EigenCompressedMatrix() = default;

    /**
     * @brief Copy constructor
     * @details The copy gets the complete row pointers even when the source is still being built by ordered
     * insertion.
     * @param rOther Matrix to copy.
     */
    EigenCompressedMatrix(const EigenCompressedMatrix& rOther) : BaseType(rOther), mAppendRow(rOther.mAppendRow)
    {
        CompleteAppend();
    }

    /**
     * @brief Move constructor
     * @details Takes over the storage and the ordered-insertion build state of rOther, and resets the build state of
     * rOther (Eigen leaves the moved-from matrix empty).
     * @param rOther Matrix whose storage and build state are taken over; it is left empty.
     */
    EigenCompressedMatrix(EigenCompressedMatrix&& rOther) noexcept : BaseType(std::move(static_cast<BaseType&>(rOther))), mAppendRow(rOther.mAppendRow)
    {
        rOther.mAppendRow = -1;
    }

    /**
     * @brief Copy assignment operator
     * @details The target gets the complete row pointers, see the copy constructor.
     * @param rOther Matrix to copy.
     * @return Reference to this matrix.
     */
    EigenCompressedMatrix& operator=(const EigenCompressedMatrix& rOther)
    {
        BaseType::operator=(rOther);
        mAppendRow = rOther.mAppendRow;
        CompleteAppend();
        return *this;
    }

    /**
     * @brief Move assignment operator
     * @details Eigen moves by swapping the storage, so the ordered-insertion build state is swapped along with it.
     * @param rOther Matrix whose storage and build state are swapped with the ones of this matrix.
     * @return Reference to this matrix.
     */
    EigenCompressedMatrix& operator=(EigenCompressedMatrix&& rOther) noexcept
    {
        // Eigen moves by swapping the storage, so the build state goes along
        BaseType::operator=(std::move(static_cast<BaseType&>(rOther)));
        std::swap(mAppendRow, rOther.mAppendRow);
        return *this;
    }

    /**
     * @brief Swap with another matrix
     * @details Swaps the storage and the ordered-insertion build state with rOther.
     * @param rOther Matrix to swap with.
     */
    void swap(EigenCompressedMatrix& rOther)
    {
        BaseType::swap(rOther);
        std::swap(mAppendRow, rOther.mAppendRow);
    }

    /**
     * @brief Constructor with the dimensions
     * @details Allocates an empty Size1 x Size2 sparse matrix (no nonzeros reserved).
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /**
     * @brief uBLAS-style (rows, cols, nnz) constructor
     * @details Allocates the compressed storage of NNZ entries up front so the CSR arrays can be written directly.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     * @param NNZ Number of entries allocated in the compressed storage.
     */
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2, const std::size_t NNZ)
        : BaseType(Size1, Size2)
    {
        this->resizeNonZeros(NNZ);
    }

    /**
     * @brief Constructor from any Eigen sparse expression
     * @details Evaluates the sparse expression rOther into this row-major matrix (implicit conversion).
     * @tparam TDerived Type of the Eigen sparse expression.
     * @param rOther Sparse expression to evaluate into this matrix.
     */
    template<class TDerived>
    EigenCompressedMatrix(const Eigen::SparseMatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /**
     * @brief Constructor from a dense Eigen expression
     * @details Explicit constructor: the nonzero entries of rExpression (taken with sparseView()) are gathered into
     * compressed storage.
     * @tparam TDerived Type of the dense Eigen expression.
     * @param rExpression Dense expression whose nonzero entries are stored.
     */
    template<class TDerived>
    explicit EigenCompressedMatrix(const Eigen::MatrixBase<TDerived>& rExpression)
        : BaseType(rExpression.derived().sparseView())
    {
        this->makeCompressed();
    }

    /**
     * @brief Constructor from a lazy zero matrix
     * @details Supports the uBLAS idiom `SparseMatrixType m = ZeroMatrix(n, n);`, which only sets the dimensions
     * (no entries are stored).
     * @param rZero Lazy zero matrix that gives the dimensions.
     */
    EigenCompressedMatrix(const zero_matrix<TDataType>& rZero) : BaseType(rZero.size1(), rZero.size2()) {}

    /**
     * @brief Constructor from a lazy identity matrix
     * @details Supports the uBLAS idiom `SparseMatrixType m = IdentityMatrix(n, n);`, which sets the dimensions and
     * stores the unit diagonal (min(size1(), size2()) entries).
     * @param rIdentity Lazy identity matrix that gives the dimensions; its unit diagonal is stored.
     */
    EigenCompressedMatrix(const identity_matrix<TDataType>& rIdentity) : BaseType(rIdentity.size1(), rIdentity.size2())
    {
        const std::size_t diagonal_size = std::min(size1(), size2());
        this->reserve(diagonal_size);
        for (std::size_t i = 0; i < diagonal_size; ++i) {
            push_back(i, i, TDataType(1));
        }
        CompleteAppend();
    }

    /**
     * @brief Conversion constructor from a compressed matrix of another scalar or index type
     * @details Explicit element-wise conversion: the stored entries of rOther are cast to TDataType.
     * @tparam TOtherDataType Scalar type of the source matrix.
     * @tparam TOtherIndexType Index type of the source matrix.
     * @param rOther Compressed matrix to convert from.
     */
    template<class TOtherDataType, class TOtherIndexType>
    explicit EigenCompressedMatrix(const EigenCompressedMatrix<TOtherDataType, TOtherIndexType>& rOther)
        : BaseType(rOther.template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment from any Eigen sparse expression
     * @details Evaluates the sparse expression rOther into this matrix and resets the ordered-insertion build state.
     * @tparam TDerived Type of the Eigen sparse expression.
     * @param rOther Sparse expression to evaluate into this matrix.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenCompressedMatrix& operator=(const Eigen::SparseMatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        mAppendRow = -1;
        return *this;
    }

    /**
     * @brief Reference to a stored entry (ublas sparse_matrix_element)
     * @details Returned by the non-const operator(), so the uBLAS idioms
     * `A(i, j) = v`, `A(i, j) += v`, `double& r = A(i, j)` and
     * `AtomicAdd(A(i, j).ref(), v)` compile unchanged.
     */
    class ElementReference
    {
    public:
        /**
         * @brief Constructor from a storage entry
         * @details Binds the proxy to rEntry, an entry of the CSR value array. The proxy refers to the entry in place,
         * so it is only valid until the next insertion into the matrix, which may move the entries.
         * @param rEntry Entry of the CSR value array the proxy refers to.
         */
        explicit ElementReference(TDataType& rEntry) : mrEntry(rEntry) {}

        /**
         * @brief Referenced storage entry
         * @details Returns a reference to the stored value, as ublas sparse_matrix_element::ref().
         * @return Reference to the stored entry.
         */
        TDataType& ref() const { return mrEntry; }

        /**
         * @brief Conversion to a reference to the stored value
         * @details Implicit conversion, so the proxy can be read or bound as a TDataType&, e.g. `double& r = A(i, j)`.
         * @return Reference to the stored entry.
         */
        operator TDataType&() const { return mrEntry; }

        /**
         * @brief Assignment of a value to the referenced entry
         * @details Writes rValue into the stored entry and returns the proxy.
         * @param rValue Value to store in the entry.
         * @return Reference to this proxy.
         */
        ElementReference& operator=(const TDataType& rValue) { mrEntry = rValue; return *this; }

        /**
         * @brief Assignment from another entry reference
         * @details Copies the value stored in the entry referenced by rOther into the entry referenced by this proxy;
         * the proxy itself is not rebound.
         * @param rOther Entry reference whose stored value is copied.
         * @return Reference to this proxy.
         */
        ElementReference& operator=(const ElementReference& rOther) { mrEntry = rOther.mrEntry; return *this; }

        /**
         * @brief Addition assignment
         * @details Adds rValue to the stored entry and returns the proxy.
         * @param rValue Value added to the entry.
         * @return Reference to this proxy.
         */
        ElementReference& operator+=(const TDataType& rValue) { mrEntry += rValue; return *this; }

        /**
         * @brief Subtraction assignment
         * @details Subtracts rValue from the stored entry and returns the proxy.
         * @param rValue Value subtracted from the entry.
         * @return Reference to this proxy.
         */
        ElementReference& operator-=(const TDataType& rValue) { mrEntry -= rValue; return *this; }

        /**
         * @brief Multiplication assignment
         * @details Multiplies the stored entry by rValue and returns the proxy.
         * @param rValue Factor the entry is multiplied by.
         * @return Reference to this proxy.
         */
        ElementReference& operator*=(const TDataType& rValue) { mrEntry *= rValue; return *this; }

        /**
         * @brief Division assignment
         * @details Divides the stored entry by rValue and returns the proxy.
         * @param rValue Divisor the entry is divided by.
         * @return Reference to this proxy.
         */
        ElementReference& operator/=(const TDataType& rValue) { mrEntry /= rValue; return *this; }

    private:
        /**
         * @brief Reference to the stored value
         * @details Entry of the CSR value array this proxy stands for.
         */
        TDataType& mrEntry;
    };

    /**
     * @brief uBLAS-style element access
     * @details As for ublas::compressed_matrix, a new entry of the last filled row (or of a later one) only shifts the
     * entries of that row, so a row-by-row construction is O(1) (amortized) per entry when the columns come in
     * increasing order; any other new entry is the O(nnz) slow path. The storage stays a packed CSR in both cases
     * (Eigen's own insertion would leave it in uncompressed mode). Returns a reference to the entry (I, J), which is
     * inserted with a zero value when it is not stored yet.
     * @param I Row index.
     * @param J Column index.
     * @return Proxy that references the entry (I, J).
     */
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

    /**
     * @brief uBLAS-style element access (const version)
     * @details Returns by value the entry (I, J), which is zero when it is not stored. The matrix is not modified.
     * @param I Row index.
     * @param J Column index.
     * @return Value of the entry (I, J), zero when it is not stored.
     */
    TDataType operator()(const std::size_t I, const std::size_t J) const { return this->coeff(I, J); }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Number of rows
     * @details uBLAS-style name of rows(), returned as std::size_t.
     * @return Number of rows.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns
     * @details uBLAS-style name of cols(), returned as std::size_t.
     * @return Number of columns.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /**
     * @brief Number of stored entries
     * @details The *filled* count (nonZeros()), as in ublas::compressed_matrix.
     * @return Number of stored entries.
     */
    std::size_t nnz() const { return static_cast<std::size_t>(this->nonZeros()); }

    /**
     * @brief Allocated storage capacity for nonzero entries
     * @details As ublas nnz_capacity(). It is the size of the compressed storage arrays (data().size()), which can
     * exceed nnz() until the row pointers declare the filled size.
     * @return Number of entries the compressed storage arrays can hold.
     */
    std::size_t nnz_capacity() const { return static_cast<std::size_t>(this->data().size()); }

    // The storage views span the *allocated* capacity (as the uBLAS storage
    // arrays do), so they can be written before the row pointers declare the
    // filled size.

    /**
     * @brief View over the CSR value array
     * @details Writable view of the stored values. It spans the allocated capacity, not only the filled entries.
     * @return Writable view over the CSR value array.
     */
    auto value_data() { return Kratos::span<TDataType>(this->valuePtr(), this->data().size()); }

    /**
     * @brief View over the CSR row-pointer array
     * @details Writable view of the row pointers (size1() + 1 entries), as ublas index1_data().
     * @return Writable view over the CSR row-pointer array.
     */
    auto index1_data() { return Kratos::span<TIndexType>(this->outerIndexPtr(), size1() + 1); }

    /**
     * @brief View over the CSR column-index array
     * @details Writable view of the column indices, as ublas index2_data(). It spans the allocated capacity, not only
     * the filled entries.
     * @return Writable view over the CSR column-index array.
     */
    auto index2_data() { return Kratos::span<TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /**
     * @brief View over the CSR value array (const version)
     * @details Read-only view of the stored values. It spans the allocated capacity, not only the filled entries.
     * @return Read-only view over the CSR value array.
     */
    auto value_data() const { return Kratos::span<const TDataType>(this->valuePtr(), this->data().size()); }

    /**
     * @brief View over the CSR row-pointer array (const version)
     * @details Read-only view of the row pointers (size1() + 1 entries).
     * @return Read-only view over the CSR row-pointer array.
     */
    auto index1_data() const { return Kratos::span<const TIndexType>(this->outerIndexPtr(), size1() + 1); }

    /**
     * @brief View over the CSR column-index array (const version)
     * @details Read-only view of the column indices. It spans the allocated capacity, not only the filled entries.
     * @return Read-only view over the CSR column-index array.
     */
    auto index2_data() const { return Kratos::span<const TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /**
     * @brief uBLAS-style filled row-pointer count
     * @details Returns size1() + 1, the number of row pointers once the matrix is assembled.
     * @return Number of row pointers, size1() + 1.
     */
    std::size_t filled1() const { return size1() + 1; }

    /**
     * @brief uBLAS-style filled nonzero count
     * @details Returns nnz().
     * @return Number of stored entries, nnz().
     */
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
     * @brief Iterator over the stored entries of a single row (ublas iterator2)
     * @details Walks the CSR arrays of one row; index1() is the row, index2() the column of the entry it
     * currently points at.
     * @tparam TIsConst Selects the read-only iterator (true), which yields const references to the values, or the
     * writable one (false).
     */
    template<bool TIsConst>
    class EntryIterator
    {
    public:
        /// Type of the pointer to the matrix. Pointer to const matrix for the read-only iterator (TIsConst), pointer to
        /// the matrix otherwise.
        using MatrixPointerType = std::conditional_t<TIsConst, const EigenCompressedMatrix*, EigenCompressedMatrix*>;

        /// Type of the reference to a value. Reference to const TDataType for the read-only iterator (TIsConst),
        /// reference to TDataType otherwise. It is the type returned by operator*.
        using reference = std::conditional_t<TIsConst, const TDataType&, TDataType&>;

        /// Scalar type of the values. Alias of TDataType.
        using value_type = TDataType;

        /**
         * @brief Constructor
         * @details Builds an iterator over the entries of the row Row of the matrix pMatrix, pointing at the position
         * Position of the CSR value and column arrays.
         * @param pMatrix Pointer to the matrix whose entries are walked.
         * @param Row Index of the row the iterator walks.
         * @param Position Position in the CSR value and column arrays the iterator points at.
         */
        EntryIterator(MatrixPointerType pMatrix, const std::size_t Row, const std::size_t Position)
            : mpMatrix(pMatrix), mRow(Row), mPosition(Position) {}

        /**
         * @brief Row of the entry currently pointed at
         * @details Returns the index of the row this iterator walks (ublas index1()).
         * @return Index of the row this iterator walks.
         */
        std::size_t index1() const { return mRow; }

        /**
         * @brief Column of the entry currently pointed at
         * @details Reads the column from the CSR column-index array at the current position (ublas index2()).
         * @return Column index of the entry currently pointed at.
         */
        std::size_t index2() const { return static_cast<std::size_t>(mpMatrix->innerIndexPtr()[mPosition]); }

        /**
         * @brief Dereference operator
         * @details Returns a reference to the value of the entry currently pointed at.
         * @return Reference to the value of the entry currently pointed at.
         */
        reference operator*() const { return mpMatrix->valuePtr()[mPosition]; }

        /**
         * @brief Pre-increment operator
         * @details Moves to the next stored entry and returns the iterator itself.
         * @return Reference to this iterator, already advanced.
         */
        EntryIterator& operator++() { ++mPosition; return *this; }

        /**
         * @brief Post-increment operator
         * @details Moves to the next stored entry and returns a copy of the iterator from before the increment.
         * @return Copy of the iterator from before the increment.
         */
        EntryIterator operator++(int) { EntryIterator copy(*this); ++mPosition; return copy; }

        /**
         * @brief Equality operator
         * @details Two iterators are equal when they point at the same position of the CSR arrays; the row is not
         * compared.
         * @param rOther Iterator to compare with.
         * @return True if both iterators point at the same position of the CSR arrays, false otherwise.
         */
        bool operator==(const EntryIterator& rOther) const { return mPosition == rOther.mPosition; }

        /**
         * @brief Inequality operator
         * @details Negation of operator==.
         * @param rOther Iterator to compare with.
         * @return True if the iterators point at different positions of the CSR arrays, false otherwise.
         */
        bool operator!=(const EntryIterator& rOther) const { return !(*this == rOther); }

    private:
        /**
         * @brief Pointer to the matrix being traversed
         * @details Pointer to const matrix for the read-only iterator.
         */
        MatrixPointerType mpMatrix;

        /**
         * @brief Row of the entries being walked
         * @details Returned by index1().
         */
        std::size_t mRow;

        /**
         * @brief Current position in the CSR arrays
         * @details Index into the CSR value/column arrays.
         */
        std::size_t mPosition;
    };

    /**
     * @brief Iterator over the rows of the matrix (ublas iterator1)
     * @details begin()/end() yield the entry iterators of the current row.
     * @tparam TIsConst Selects the read-only iterator (true), whose entry iterators yield const references to the
     * values, or the writable one (false).
     */
    template<bool TIsConst>
    class RowIterator
    {
    public:
        /// Type of the pointer to the matrix. Pointer to const matrix for the read-only iterator (TIsConst), pointer to
        /// the matrix otherwise.
        using MatrixPointerType = std::conditional_t<TIsConst, const EigenCompressedMatrix*, EigenCompressedMatrix*>;

        /// Type of the iterator over the entries of a row. EntryIterator with the same constness as this iterator,
        /// returned by begin() and end().
        using EntryIteratorType = EntryIterator<TIsConst>;

        /**
         * @brief Constructor
         * @details Builds an iterator pointing at the row Row of the matrix pMatrix.
         * @param pMatrix Pointer to the matrix whose rows are walked.
         * @param Row Index of the row the iterator points at.
         */
        RowIterator(MatrixPointerType pMatrix, const std::size_t Row)
            : mpMatrix(pMatrix), mRow(Row) {}

        /**
         * @brief Index of the row currently pointed at
         * @details Returns the row index (ublas index1()).
         * @return Index of the row currently pointed at.
         */
        std::size_t index1() const { return mRow; }

        /**
         * @brief First stored entry of this row
         * @details Returns an entry iterator at the position given by the row pointer of this row.
         * @return Entry iterator at the first stored entry of this row.
         */
        EntryIteratorType begin() const
        {
            return EntryIteratorType(mpMatrix, mRow, static_cast<std::size_t>(mpMatrix->outerIndexPtr()[mRow]));
        }

        /**
         * @brief Past-the-last stored entry of this row
         * @details Returns an entry iterator at the position given by the row pointer of the next row.
         * @return Entry iterator past the last stored entry of this row.
         */
        EntryIteratorType end() const
        {
            return EntryIteratorType(mpMatrix, mRow, static_cast<std::size_t>(mpMatrix->outerIndexPtr()[mRow + 1]));
        }

        /**
         * @brief Pre-increment operator
         * @details Moves to the next row and returns the iterator itself.
         * @return Reference to this iterator, already advanced.
         */
        RowIterator& operator++() { ++mRow; return *this; }

        /**
         * @brief Post-increment operator
         * @details Moves to the next row and returns a copy of the iterator from before the increment.
         * @return Copy of the iterator from before the increment.
         */
        RowIterator operator++(int) { RowIterator copy(*this); ++mRow; return copy; }

        /**
         * @brief Equality operator
         * @details Two iterators are equal when they point at the same row.
         * @param rOther Iterator to compare with.
         * @return True if both iterators point at the same row, false otherwise.
         */
        bool operator==(const RowIterator& rOther) const { return mRow == rOther.mRow; }

        /**
         * @brief Inequality operator
         * @details Negation of operator==.
         * @param rOther Iterator to compare with.
         * @return True if the iterators point at different rows, false otherwise.
         */
        bool operator!=(const RowIterator& rOther) const { return !(*this == rOther); }

    private:
        /**
         * @brief Pointer to the matrix being traversed
         * @details Pointer to const matrix for the read-only iterator.
         */
        MatrixPointerType mpMatrix;

        /**
         * @brief Row currently pointed at
         * @details Returned by index1().
         */
        std::size_t mRow;
    };

    /// Writable iterator over the rows. `RowIterator<false>`, the uBLAS iterator1.
    using iterator1 = RowIterator<false>;

    /// Read-only iterator over the rows. `RowIterator<true>`, the uBLAS const_iterator1.
    using const_iterator1 = RowIterator<true>;

    /// Writable iterator over the entries of a row. `EntryIterator<false>`, the uBLAS iterator2.
    using iterator2 = EntryIterator<false>;

    /// Read-only iterator over the entries of a row. `EntryIterator<true>`, the uBLAS const_iterator2.
    using const_iterator2 = EntryIterator<true>;

    /**
     * @brief First row
     * @details The storage is packed first (complete_index1_data()), so the CSR arrays the iterators walk are the
     * compressed ones; the packing is a no-op when the storage already is compressed. Returns a writable row
     * iterator at the row 0.
     * @return Writable row iterator at the row 0.
     */
    iterator1 begin1()
    {
        complete_index1_data();
        return iterator1(this, 0);
    }

    /**
     * @brief Past-the-last row
     * @details Returns a writable row iterator at the row size1().
     * @return Writable row iterator at the row size1().
     */
    iterator1 end1() { return iterator1(this, size1()); }

    /**
     * @brief First row (const version)
     * @details Returns a read-only row iterator at the row 0. Being const it does not pack the storage: in debug
     * builds an error is raised if the storage is not compressed, so call complete_index1_data() first.
     * @return Read-only row iterator at the row 0.
     */
    const_iterator1 begin1() const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(this->isCompressed())
            << "Iterating a compressed matrix whose storage is not packed; call complete_index1_data() first." << std::endl;
        return const_iterator1(this, 0);
    }

    /**
     * @brief Past-the-last row (const version)
     * @details Returns a read-only row iterator at the row size1().
     * @return Read-only row iterator at the row size1().
     */
    const_iterator1 end1() const { return const_iterator1(this, size1()); }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief uBLAS-style finalization after writing the CSR arrays by hand
     * @details The storage was already sized and the filled count is implied by the written row pointers, so this
     * only validates consistency: it completes the pending row pointers and, in debug builds, raises an error if
     * FilledSize1 is not size1() + 1, if FilledNNZ differs from nnz() or if FilledNNZ exceeds the allocated capacity.
     * @param FilledSize1 Number of row pointers declared as filled; it must be size1() + 1.
     * @param FilledNNZ Number of entries declared as filled; it must equal nnz() and not exceed nnz_capacity().
     */
    void set_filled(const std::size_t FilledSize1, const std::size_t FilledNNZ)
    {
        CompleteAppend();
        KRATOS_DEBUG_ERROR_IF(FilledSize1 != size1() + 1 || FilledNNZ != nnz() || FilledNNZ > static_cast<std::size_t>(this->data().size()))
            << "set_filled(" << FilledSize1 << ", " << FilledNNZ
            << ") is inconsistent with the written compressed storage ("
            << size1() + 1 << ", " << nnz() << " of " << this->data().size()
            << " allocated)." << std::endl;
    }

    /**
     * @brief uBLAS-style resize
     * @details As in ublas::compressed_matrix the default preserves, while resize(m, n, false) discards values and
     * structure. The pending row pointers are completed first, then Preserve selects Eigen's conservativeResize
     * (true) or resize (false).
     * @param NewSize1 New number of rows.
     * @param NewSize2 New number of columns.
     * @param Preserve If true (default) the existing values and structure are kept (Eigen's conservativeResize), if
     * false they are discarded (Eigen's resize).
     */
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

    /**
     * @brief uBLAS-style reservation of nonzero storage capacity
     * @details Preserve=true keeps the existing structure and calls Eigen's reserve. Preserve=false only resizes the
     * compressed storage to NNZ entries (resizeNonZeros) and does not clear the entries already stored. The pending row
     * pointers are completed first.
     * @param NNZ Number of nonzero entries to reserve (the size of the compressed storage when Preserve is false).
     * @param Preserve If true (default) the existing structure is kept and Eigen's reserve is called, if false only the
     * compressed storage is resized to NNZ entries.
     */
    void reserve(const std::size_t NNZ, const bool Preserve = true)
    {
        CompleteAppend();
        if (Preserve) {
            BaseType::reserve(static_cast<typename BaseType::Index>(NNZ));
        } else {
            this->resizeNonZeros(static_cast<typename BaseType::Index>(NNZ));
        }
    }

    /**
     * @brief uBLAS-style clear
     * @details Removes all stored entries (the dimensions are kept) and resets the ordered-insertion build state.
     */
    void clear()
    {
        mAppendRow = -1;
        this->setZero();
    }

    /**
     * @brief uBLAS-style ordered insertion
     * @details As ublas::compressed_matrix::push_back(): entries given row by row with increasing columns are
     * appended in O(1) (amortized). Sets the entry (I, J) to Value through operator().
     * @param I Row index.
     * @param J Column index.
     * @param Value Value assigned to the entry (I, J).
     */
    void push_back(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        (*this)(I, J) = Value;
    }

    /**
     * @brief uBLAS-style element insertion
     * @details As ublas::compressed_matrix::insert_element. Sets the entry (I, J) to Value through operator().
     * @param I Row index.
     * @param J Column index.
     * @param Value Value assigned to the entry (I, J).
     */
    void insert_element(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        (*this)(I, J) = Value;
    }

    /**
     * @brief uBLAS-style finalization
     * @details As ublas::compressed_matrix::complete_index1_data(): sets the row pointers of the trailing rows left
     * empty by the ordered insertion, and packs the storage after Eigen-native insertions (insert(), reserve()) left
     * it uncompressed.
     */
    void complete_index1_data()
    {
        CompleteAppend();
        this->makeCompressed();
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /**
     * @brief Build state of the ordered insertion
     * @details Row of the last entry appended by the ordered insertion while the row pointers of the rows after it
     * are not yet written (-1 otherwise). As for ublas::compressed_matrix::filled1_, the row pointers past it are
     * completed by the next non-ordered access or complete_index1_data().
     */
    std::ptrdiff_t mAppendRow = -1;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Number of entries written so far
     * @details The packed filled count. While an ordered insertion is pending it is read from the row pointer that
     * follows the last appended row, otherwise from the last row pointer.
     * @return Number of entries written so far in the compressed storage.
     */
    std::size_t FilledEntries() const
    {
        const TIndexType* p_outer = this->outerIndexPtr();
        return static_cast<std::size_t>(mAppendRow >= 0 ? p_outer[mAppendRow + 1] : p_outer[size1()]);
    }

    /**
     * @brief Check whether row I is a tail row
     * @details True if row I is the last row with stored entries or comes after it, so an entry of it can be inserted
     * by shifting that row only. Always false when the storage is not compressed or when (I, J) is out of range.
     * @param I Row index.
     * @param J Column index, only used to check that the entry (I, J) is in range.
     * @return True if an entry of row I can be inserted by shifting that row only, false otherwise.
     */
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

    /**
     * @brief Entry (I, J) of the tail row I
     * @details IsTailRow(I, J) must hold. The entry is inserted with a zero value if not stored yet. The columns of
     * the row are kept sorted, shifting the entries of that row only: an ordered insertion is an O(1) (amortized)
     * append, as for ublas::compressed_matrix. Returns a reference to the stored value.
     * @param I Row index, which must be a tail row (see IsTailRow()).
     * @param J Column index.
     * @return Reference to the stored value of the entry (I, J).
     */
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

    /**
     * @brief Complete the row pointers pending after the ordered insertion
     * @details Writes the row pointers of the rows after the last appended one (all equal to the filled count) and
     * resets the build state. It does nothing when no row pointer is pending.
     */
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
