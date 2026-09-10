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
#include <istream>
#include <ostream>
#include <type_traits>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <boost/numeric/ublas/expression_types.hpp>

// Project includes
#include "includes/exception.h"

namespace Kratos {

/**
 * @brief Index type used for the Eigen sparse matrices.
 * @details Eigen requires a signed StorageIndex. The default is a 32-bit int:
 * sparse kernels are memory-bound (SpMV traffic drops from 16 to 12 bytes per
 * nonzero versus 64-bit indices, a structural advantage the std::size_t-indexed
 * uBLAS backend cannot match), and 2^31 - 1 nonzeros per matrix (~16 GB of
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
static_assert(std::is_signed_v<KratosEigenIndexType>,
              "Eigen requires a signed sparse StorageIndex.");

namespace Internals {

/**
 * @class EigenArrayProxy
 * @brief Iterable view over a raw backend array.
 * @details Gives Eigen's CSR storage arrays the same begin()/end()/operator[]
 * surface that the uBLAS unbounded_array storage exposes, so code written
 * against compressed_matrix::value_data()/index1_data()/index2_data() works
 * unchanged on the Eigen wrapper types.
 * @tparam T The (possibly const-qualified) element type of the viewed array.
 */
template<class T>
class EigenArrayProxy
{
public:
    using value_type = std::remove_const_t<T>; /// Element type, with any const-qualification of T stripped.
    using iterator = T*; /// Iterator type (a raw pointer, const-qualified iff T is).
    using const_iterator = const T*; /// Const iterator type.

    /**
     * @brief Constructor wrapping a raw pointer and its length.
     * @param pData Pointer to the first element of the backend array.
     * @param Size Number of elements the array holds.
     */
    EigenArrayProxy(T* pData, const std::size_t Size) : mpData(pData), mSize(Size) {}

    /**
     * @brief Returns an iterator to the first element.
     * @return Pointer to the first element.
     */
    T* begin() { return mpData; }

    /**
     * @brief Returns an iterator past the last element.
     * @return Pointer past the last element.
     */
    T* end() { return mpData + mSize; }

    /**
     * @brief Returns a const iterator to the first element.
     * @return Const pointer to the first element.
     */
    const T* begin() const { return mpData; }

    /**
     * @brief Returns a const iterator past the last element.
     * @return Const pointer past the last element.
     */
    const T* end() const { return mpData + mSize; }

    /**
     * @brief Element access by index.
     * @param Index Zero-based index of the requested element.
     * @return Reference to the element at Index.
     */
    T& operator[](const std::size_t Index) { return mpData[Index]; }

    /**
     * @brief Element access by index (const version).
     * @param Index Zero-based index of the requested element.
     * @return Const reference to the element at Index.
     */
    const T& operator[](const std::size_t Index) const { return mpData[Index]; }

    /**
     * @brief Number of elements viewed by this proxy.
     * @return The array size.
     */
    std::size_t size() const { return mSize; }

private:
    T* mpData; /// Pointer to the first element of the viewed backend array.
    std::size_t mSize; /// Number of elements viewed.
};

} // namespace Internals

/**
 * @class EigenMatrix
 * @brief Dense, dynamically sized, row-major Eigen matrix with the uBLAS member surface.
 * @details Row-major storage is chosen on purpose: it matches the storage
 * order of the Kratos (uBLAS) dense matrices, so raw-buffer interoperability
 * is preserved. The added members (size1/size2 and the 3-argument resize)
 * mirror boost::numeric::ublas::matrix so generic Kratos code compiles
 * unchanged against this type.
 * @tparam TDataType The scalar type stored in the matrix (e.g. double).
 */
template<class TDataType>
class EigenMatrix : public Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
public:
    using BaseType = Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; /// The wrapped Eigen dense matrix type.
    using value_type = TDataType; /// Scalar type stored in the matrix.
    using size_type = std::size_t; /// uBLAS-style unsigned size type.

    // Minimal uBLAS container-trait surface, mirroring the EigenVector one:
    // the generic boost::numeric::ublas::project() overloads instantiate
    // matrix_range<M>/matrix_slice<M> while they are *considered* during
    // overload resolution, and that class instantiation needs these member
    // typedefs to be well-formed. The proxies are never actually used on this
    // type - the Kratos Eigen overloads win the resolution - so the closure
    // aliases are plain references and the (never dereferenced through the
    // uBLAS iterator protocol) iterator aliases are plain pointers.
    using difference_type = std::ptrdiff_t; /// uBLAS-style signed difference type.
    using reference = TDataType&; /// uBLAS-style mutable element reference.
    using const_reference = const TDataType&; /// uBLAS-style const element reference.
    using closure_type = EigenMatrix<TDataType>&; /// uBLAS-style closure (proxy) type; never dereferenced through the uBLAS iterator protocol, so a plain reference suffices.
    using const_closure_type = const EigenMatrix<TDataType>&; /// uBLAS-style const closure (proxy) type.
    using storage_category = boost::numeric::ublas::dense_tag; /// uBLAS storage-layout trait: dense storage.
    using orientation_category = boost::numeric::ublas::row_major_tag; /// uBLAS orientation trait: row-major storage.
    using iterator1 = TDataType*; /// uBLAS-style row (first-dimension) iterator; a plain pointer, never dereferenced through the uBLAS iterator protocol.
    using const_iterator1 = const TDataType*; /// uBLAS-style const row (first-dimension) iterator.
    using iterator2 = TDataType*; /// uBLAS-style column (second-dimension) iterator; a plain pointer, never dereferenced through the uBLAS iterator protocol.
    using const_iterator2 = const TDataType*; /// uBLAS-style const column (second-dimension) iterator.

    /**
     * @brief Default constructor. Creates an empty (0x0) matrix.
     */
    EigenMatrix() = default;

    /**
     * @brief Constructor allocating an uninitialized Size1 x Size2 matrix.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    EigenMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /**
     * @brief uBLAS-style (rows, cols, value) fill constructor.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     * @param Value Value assigned to every entry.
     */
    EigenMatrix(const std::size_t Size1, const std::size_t Size2, const TDataType Value)
        : BaseType(BaseType::Constant(Size1, Size2, Value)) {}

    /**
     * @brief Construction from any Eigen expression.
     * @details The scalar cast makes the cross-precision conversions of ublas
     * (double <-> float grids in the p-multigrid hierarchy) work; it is the
     * identity for matching scalars.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     */
    template<class TDerived>
    EigenMatrix(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    /**
     * @brief Assignment from any Eigen expression (scalar-casting, see the constructor).
     * @details The assignment evaluates through a temporary (.eval()), reproducing
     * the alias-safe uBLAS operator= semantics (e.g. M = trans(M), which Eigen
     * would otherwise reject at runtime); the temporary-free fast path remains
     * noalias(M) = expr / assign(), as in uBLAS.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenMatrix& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>().eval());
        return *this;
    }

    /**
     * @brief Construction from a dense uBLAS matrix expression.
     * @details Covers ZeroMatrix, IdentityMatrix, prod(...) results, etc.,
     * mirroring ublas::matrix's own converting constructor so uBLAS idioms
     * keep compiling on this type.
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     */
    template<class TExpression>
    EigenMatrix(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
        : BaseType(rExpression().size1(), rExpression().size2())
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                (*this)(i, j) = r_e(i, j);
    }

    /**
     * @brief Assignment from a dense uBLAS matrix expression (resizing, as ublas does).
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     * @return Reference to this matrix.
     */
    template<class TExpression>
    EigenMatrix& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        BaseType::resize(r_e.size1(), r_e.size2());
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                (*this)(i, j) = r_e(i, j);
        return *this;
    }

    /**
     * @brief uBLAS-style assign: element-wise assignment without a protective
     * temporary (ublas::matrix::assign), the operation the generic linear
     * algebra spaces call for their Copy(rX, rY).
     * @details uBLAS requires the target to be sized already; this resizes it
     * when it is not, a superset of that contract and what the Eigen space
     * does in its own Copy.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenMatrix& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS-style assign from a dense uBLAS matrix expression.
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     * @return Reference to this matrix.
     */
    template<class TExpression>
    EigenMatrix& assign(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        return operator=(rExpression);
    }

    /**
     * @brief Number of rows.
     * @return Row count.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns.
     * @return Column count.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /**
     * @brief uBLAS-style resize; note that, as in ublas::matrix, the default preserves the values.
     * @param NewSize1 New number of rows.
     * @param NewSize2 New number of columns.
     * @param Preserve If true (default) existing values are preserved; if false they are discarded.
     */
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize1, NewSize2);
        else BaseType::resize(NewSize1, NewSize2);
    }

    /**
     * @brief uBLAS-style clear: zero all entries (the size is kept).
     */
    void clear()
    {
        this->setZero();
    }
};

/**
 * @class EigenVector
 * @brief Dense, dynamically sized Eigen column vector with the uBLAS member surface.
 * @details Adds the (size), (size, value) constructors and the 2-argument
 * resize of boost::numeric::ublas::vector; operator[], operator() and size()
 * are already provided by Eigen with compatible semantics.
 * @tparam TDataType The scalar type stored in the vector (e.g. double).
 */
template<class TDataType>
class EigenVector : public Eigen::Matrix<TDataType, Eigen::Dynamic, 1>
{
public:
    using BaseType = Eigen::Matrix<TDataType, Eigen::Dynamic, 1>; /// The wrapped Eigen dense column-vector type.
    using value_type = TDataType; /// Scalar type stored in the vector.
    using size_type = std::size_t; /// uBLAS-style unsigned size type.

    // Minimal uBLAS container-trait surface. The generic
    // boost::numeric::ublas::project()/subrange() overloads instantiate
    // vector_range<V>/vector_slice<V> while they are *considered* during
    // overload resolution (their second parameter names a nested type of the
    // proxy), and that class instantiation needs these member typedefs to be
    // well-formed. The proxies are never actually used on this type - the
    // Kratos Eigen overloads win the resolution - so reference-typed closure
    // aliases are sufficient.
    using difference_type = std::ptrdiff_t; /// uBLAS-style signed difference type.
    using reference = TDataType&; /// uBLAS-style mutable element reference.
    using const_reference = const TDataType&; /// uBLAS-style const element reference.
    // iterator/const_iterator are the inherited Eigen STL iterator typedefs (TDataType*/const TDataType*).
    using closure_type = EigenVector<TDataType>&; /// uBLAS-style closure (proxy) type.
    using const_closure_type = const EigenVector<TDataType>&; /// uBLAS-style const closure (proxy) type.
    using storage_category = boost::numeric::ublas::dense_tag; /// uBLAS storage-layout trait: dense storage.

    /**
     * @brief Default constructor. Creates an empty (0-length) vector.
     */
    EigenVector() = default;

    /**
     * @brief Constructor allocating an uninitialized vector of the given size.
     * @param Size Number of entries.
     */
    explicit EigenVector(const std::size_t Size) : BaseType(Size) {}

    /**
     * @brief uBLAS-style (size, value) fill constructor.
     * @param Size Number of entries.
     * @param Value Value assigned to every entry.
     */
    EigenVector(const std::size_t Size, const TDataType Value)
        : BaseType(BaseType::Constant(Size, Value)) {}

    /**
     * @brief Construction from any Eigen expression.
     * @details The scalar cast makes the cross-precision conversions of ublas
     * (double <-> float grids in the p-multigrid hierarchy) work; it is the
     * identity for matching scalars.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     */
    template<class TDerived>
    EigenVector(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    /**
     * @brief Assignment from any Eigen expression (scalar-casting, see the constructor).
     * @details The assignment evaluates through a temporary (.eval()), reproducing
     * the alias-safe uBLAS operator= semantics; the temporary-free fast path
     * remains noalias(v) = expr / assign(), as in uBLAS.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenVector& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>().eval());
        return *this;
    }

    /**
     * @brief Construction from a dense uBLAS vector expression.
     * @details Covers ZeroVector, prod(...) results, etc., mirroring
     * ublas::vector's own converting constructor so uBLAS idioms keep
     * compiling on this type.
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     */
    template<class TExpression>
    EigenVector(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
        : BaseType(rExpression().size())
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) (*this)[i] = r_e(i);
    }

    /**
     * @brief Assignment from a dense uBLAS vector expression (resizing, as ublas does).
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     * @return Reference to this vector.
     */
    template<class TExpression>
    EigenVector& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        BaseType::resize(r_e.size());
        for (std::size_t i = 0; i < r_e.size(); ++i) (*this)[i] = r_e(i);
        return *this;
    }

    /**
     * @brief uBLAS-style assign: element-wise assignment without a protective
     * temporary (ublas::vector::assign), the operation the generic linear
     * algebra spaces call for their Copy(rX, rY).
     * @details uBLAS requires the target to be sized already; this resizes it
     * when it is not, a superset of that contract and what the Eigen space
     * does in its own Copy.
     * @tparam TDerived Concrete type of the Eigen expression.
     * @param rOther Source Eigen expression.
     * @return Reference to this vector.
     */
    template<class TDerived>
    EigenVector& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /**
     * @brief uBLAS-style assign from a dense uBLAS vector expression.
     * @tparam TExpression Concrete type of the uBLAS vector expression.
     * @param rExpression Source uBLAS vector expression.
     * @return Reference to this vector.
     */
    template<class TExpression>
    EigenVector& assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return operator=(rExpression);
    }

    /**
     * @brief uBLAS-style resize; note that, as in ublas::vector, the default preserves the values.
     * @param NewSize New number of entries.
     * @param Preserve If true (default) existing values are preserved; if false they are discarded.
     */
    void resize(const std::size_t NewSize, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize);
        else BaseType::resize(NewSize);
    }

    /**
     * @brief uBLAS-style clear: zero all entries (the size is kept).
     */
    void clear()
    {
        this->setZero();
    }

    /**
     * @brief uBLAS-style size: unsigned, so the pervasive comparisons/loops
     * against std::size_t counters stay warning-free (Eigen's own size() is signed).
     * @return Number of entries.
     */
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }

    /**
     * @brief uBLAS-style emptiness check.
     * @return True if the vector has no entries.
     */
    bool empty() const { return BaseType::size() == 0; }
};

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
    using BaseType = Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>; /// The wrapped Eigen row-major sparse matrix type.
    using value_type = TDataType; /// Scalar type stored in the matrix.
    using size_type = std::size_t; /// uBLAS-style unsigned size type.
    // uBLAS-style storage-array typedefs (generic code uses e.g.
    // typename TMatrix::index_array_type::value_type to name the index type)
    using index_array_type = Internals::EigenArrayProxy<TIndexType>; /// Type of the index1_data()/index2_data() CSR index-array proxies.
    using value_array_type = Internals::EigenArrayProxy<TDataType>; /// Type of the value_data() CSR value-array proxy.

    /**
     * @brief Default constructor. Creates an empty (0x0) matrix.
     */
    EigenCompressedMatrix() = default;

    /**
     * @brief Constructor allocating an empty Size1 x Size2 sparse matrix (no nonzeros reserved).
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     */
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /**
     * @brief uBLAS-style (rows, cols, nnz) constructor.
     * @details Allocates the compressed storage up front so the CSR arrays
     * can be written directly.
     * @param Size1 Number of rows.
     * @param Size2 Number of columns.
     * @param NNZ Number of nonzero entries to reserve storage for.
     */
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2, const std::size_t NNZ)
        : BaseType(Size1, Size2)
    {
        this->resizeNonZeros(NNZ);
    }

    /**
     * @brief Construction from any Eigen sparse expression.
     * @tparam TDerived Concrete type of the Eigen sparse expression.
     * @param rOther Source Eigen sparse expression.
     */
    template<class TDerived>
    EigenCompressedMatrix(const Eigen::SparseMatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /**
     * @brief Construction from a dense Eigen expression.
     * @details The nonzero entries are gathered into compressed storage (the
     * dense counterpart of the sparse constructor above), mirroring
     * ublas::compressed_matrix's converting constructor for the Eigen-backed
     * dense types.
     * @tparam TDerived Concrete type of the Eigen dense expression.
     * @param rExpression Source Eigen dense expression.
     */
    template<class TDerived>
    explicit EigenCompressedMatrix(const Eigen::MatrixBase<TDerived>& rExpression)
        : BaseType(rExpression.derived().sparseView())
    {
        this->makeCompressed();
    }

    /**
     * @brief Construction from a dense uBLAS matrix expression.
     * @details Mirrors ublas::compressed_matrix's own constructor: the
     * nonzero entries are gathered into compressed storage. Insertion happens
     * in row-major order (this matrix's storage order), so it is a single
     * linear fill. Not explicit: EigenMatrix's equivalent constructor isn't
     * either, and the ublas idiom `SparseMatrixType m = ZeroMatrix(n, n);`
     * relies on the implicit conversion (as does its uBLAS-backend counterpart).
     * @tparam TExpression Concrete type of the uBLAS matrix expression.
     * @param rExpression Source uBLAS matrix expression.
     */
    template<class TExpression>
    EigenCompressedMatrix(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
        : BaseType(rExpression().size1(), rExpression().size2())
    {
        const auto& r_dense = rExpression();
        const std::size_t size1 = r_dense.size1();
        const std::size_t size2 = r_dense.size2();

        Eigen::Matrix<TIndexType, Eigen::Dynamic, 1> nnz_per_row = Eigen::Matrix<TIndexType, Eigen::Dynamic, 1>::Zero(size1);
        for (std::size_t i = 0; i < size1; ++i)
            for (std::size_t j = 0; j < size2; ++j)
                if (r_dense(i, j) != TDataType{}) ++nnz_per_row[i];

        this->reserve(nnz_per_row);
        for (std::size_t i = 0; i < size1; ++i)
            for (std::size_t j = 0; j < size2; ++j)
                if (r_dense(i, j) != TDataType{}) this->insert(i, j) = r_dense(i, j);
        this->makeCompressed();
    }

    /**
     * @brief Assignment from any Eigen sparse expression.
     * @tparam TDerived Concrete type of the Eigen sparse expression.
     * @param rOther Source Eigen sparse expression.
     * @return Reference to this matrix.
     */
    template<class TDerived>
    EigenCompressedMatrix& operator=(const Eigen::SparseMatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /**
     * @brief Number of rows.
     * @return Row count.
     */
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /**
     * @brief Number of columns.
     * @return Column count.
     */
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /**
     * @brief Number of stored entries.
     * @details As in ublas::compressed_matrix this is the *filled* count (in
     * Eigen's compressed mode it is derived from the row pointers), which is
     * 0 until the CSR arrays are written; the allocated capacity is exposed
     * through the size of the value/index proxies.
     * @return Number of filled (stored) nonzero entries.
     */
    std::size_t nnz() const { return static_cast<std::size_t>(this->nonZeros()); }

    /**
     * @brief Allocated storage capacity for nonzero entries.
     * @details As in ublas::compressed_matrix::nnz_capacity(), this is the
     * size of the allocated value_data()/index2_data() storage, which may
     * exceed nnz() until the compressed storage is fully written.
     * @return Allocated nonzero capacity.
     */
    std::size_t nnz_capacity() const { return static_cast<std::size_t>(this->data().size()); }

    // The storage proxies span the *allocated* capacity (as the uBLAS storage
    // arrays do), so they can be written before the row pointers declare the
    // filled size.

    /**
     * @brief Proxy over the CSR value array.
     * @return Iterable view spanning the allocated value storage.
     */
    auto value_data() { return Internals::EigenArrayProxy<TDataType>(this->valuePtr(), this->data().size()); }

    /**
     * @brief Proxy over the CSR row-pointer array.
     * @return Iterable view spanning the size1() + 1 row pointers.
     */
    auto index1_data() { return Internals::EigenArrayProxy<TIndexType>(this->outerIndexPtr(), size1() + 1); }

    /**
     * @brief Proxy over the CSR column-index array.
     * @return Iterable view spanning the allocated column-index storage.
     */
    auto index2_data() { return Internals::EigenArrayProxy<TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /**
     * @brief Proxy over the CSR value array (const version).
     * @return Const iterable view spanning the allocated value storage.
     */
    auto value_data() const { return Internals::EigenArrayProxy<const TDataType>(this->valuePtr(), this->data().size()); }

    /**
     * @brief Proxy over the CSR row-pointer array (const version).
     * @return Const iterable view spanning the size1() + 1 row pointers.
     */
    auto index1_data() const { return Internals::EigenArrayProxy<const TIndexType>(this->outerIndexPtr(), size1() + 1); }

    /**
     * @brief Proxy over the CSR column-index array (const version).
     * @return Const iterable view spanning the allocated column-index storage.
     */
    auto index2_data() const { return Internals::EigenArrayProxy<const TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /**
     * @brief uBLAS-style filled row-pointer count.
     * @details Mirrors ublas::compressed_matrix::filled1(): the number of
     * written entries in index1_data(), i.e. size1() + 1 once fully assembled.
     * @return Number of filled row-pointer entries.
     */
    std::size_t filled1() const { return size1() + 1; }

    /**
     * @brief uBLAS-style filled nonzero count.
     * @details Mirrors ublas::compressed_matrix::filled2(): the number of
     * written entries in index2_data()/value_data(), i.e. nnz().
     * @return Number of filled nonzero entries.
     */
    std::size_t filled2() const { return nnz(); }

    /**
     * @brief uBLAS-style finalization after writing the CSR arrays by hand.
     * @details The storage was already sized by the (rows, cols, nnz)
     * constructor and the filled count is implied by the written row
     * pointers, so this only validates consistency.
     * @param FilledSize1 Expected number of written row pointers (size1() + 1).
     * @param FilledNNZ Expected number of written nonzero entries.
     */
    void set_filled(const std::size_t FilledSize1, const std::size_t FilledNNZ)
    {
        KRATOS_DEBUG_ERROR_IF(FilledSize1 != size1() + 1 || FilledNNZ != nnz() || FilledNNZ > static_cast<std::size_t>(this->data().size()))
            << "set_filled(" << FilledSize1 << ", " << FilledNNZ
            << ") is inconsistent with the written compressed storage ("
            << size1() + 1 << ", " << nnz() << " of " << this->data().size()
            << " allocated)." << std::endl;
    }

    /**
     * @brief uBLAS-style resize; as in ublas::compressed_matrix the default preserves,
     * while resize(m, n, false) discards both values and structure.
     * @param NewSize1 New number of rows.
     * @param NewSize2 New number of columns.
     * @param Preserve If true (default) existing values and structure are preserved; if false they are discarded.
     */
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize1, NewSize2);
        else BaseType::resize(NewSize1, NewSize2);
    }

    // Brings in BaseType's reserve() overloads (single count, or per-inner-vector
    // counts) alongside the uBLAS-style two-argument one added below.
    using BaseType::reserve;

    /**
     * @brief uBLAS-style reservation of nonzero storage capacity.
     * @details As in ublas::compressed_matrix::reserve, Preserve=false
     * discards existing values and structure while allocating room for NNZ
     * nonzeros (mirroring the (rows, cols, nnz) constructor); Preserve=true
     * keeps the existing structure and behaves as the base class's own reserve.
     * @param NNZ Number of nonzero entries to allocate storage for.
     * @param Preserve If true (default) existing values and structure are preserved; if false they are discarded.
     */
    void reserve(const std::size_t NNZ, const bool Preserve = true)
    {
        if (Preserve) BaseType::reserve(static_cast<typename BaseType::Index>(NNZ));
        else this->resizeNonZeros(static_cast<typename BaseType::Index>(NNZ));
    }

    /**
     * @brief uBLAS-style clear: removes all stored entries (the dimensions are kept).
     */
    void clear()
    {
        this->setZero();
    }

    /**
     * @brief uBLAS-style element insertion.
     * @details As for operator() on a missing entry this is an O(nnz-in-row)
     * slow path meant for tests and small setup code, not for assembly
     * (write the CSR arrays directly instead).
     * @param I Row index.
     * @param J Column index.
     * @param Value Value to store at (I, J).
     */
    void push_back(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        this->coeffRef(I, J) = Value;
    }

    /**
     * @brief uBLAS-style finalization after element-wise insertion, as
     * ublas::compressed_matrix::complete_index1_data().
     * @details Element insertion (push_back, insert_element, operator() on a
     * missing entry) leaves the Eigen storage uncompressed, where
     * index1_data()/index2_data()/value_data() are not the packed CSR arrays.
     * Call this once the pattern is built and before reading those arrays.
     */
    void complete_index1_data()
    {
        this->makeCompressed();
    }

    /**
     * @brief uBLAS-style element insertion, as ublas::compressed_matrix::insert_element.
     * @details Unlike push_back it does not require appending in row-major
     * order, at the same O(nnz-in-row) slow-path cost as operator() on a
     * missing entry.
     * @param I Row index.
     * @param J Column index.
     * @param Value Value to store at (I, J).
     */
    void insert_element(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        this->coeffRef(I, J) = Value;
    }

    /**
     * @brief uBLAS-style element access.
     * @details Inserting a new entry through the non-const version is the
     * same slow path as for ublas::compressed_matrix.
     * @param I Row index.
     * @param J Column index.
     * @return Reference to the entry at (I, J), inserting it if missing.
     */
    TDataType& operator()(const std::size_t I, const std::size_t J) { return this->coeffRef(I, J); }

    /**
     * @brief uBLAS-style element access (const version).
     * @param I Row index.
     * @param J Column index.
     * @return Value of the entry at (I, J), or zero if not stored.
     */
    TDataType operator()(const std::size_t I, const std::size_t J) const { return this->coeff(I, J); }
};

/**@name uBLAS-format stream operators */
/**@{*/
// The dynamic Eigen-backed types print and parse in the boost::numeric::ublas
// text format ("[n](v0,v1,...)" and "[n1,n2]((a00,a01),(a10,a11))"): the
// format is part of the established Kratos IO surface (mdpa Matrix/Vector
// values, json checks, restart files, printed expectations in tests), so it
// must not change with the backend. These exact-type overloads win over
// Eigen's generic MatrixBase operator<<.

/**
 * @brief Writes an EigenVector in the uBLAS text format "[n](v0,v1,...)".
 * @tparam TDataType The vector's scalar type.
 * @param rOStream Output stream to write to.
 * @param rV Vector to serialize.
 * @return Reference to rOStream.
 */
template<class TDataType>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenVector<TDataType>& rV)
{
    rOStream << '[' << rV.size() << "](";
    for (std::size_t i = 0; i < static_cast<std::size_t>(rV.size()); ++i) {
        if (i != 0) rOStream << ',';
        rOStream << rV[i];
    }
    rOStream << ')';
    return rOStream;
}

/**
 * @brief Reads an EigenVector from the uBLAS text format "[n](v0,v1,...)".
 * @details Parses the format written by the matching operator<<. On any
 * mismatch the failbit is set, as boost's own operator>> does.
 * @tparam TDataType The vector's scalar type.
 * @param rIStream Input stream to read from.
 * @param rV Vector to populate; resized to match the parsed size.
 * @return Reference to rIStream.
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
 * @brief Writes an EigenMatrix in the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))".
 * @tparam TDataType The matrix's scalar type.
 * @param rOStream Output stream to write to.
 * @param rM Matrix to serialize.
 * @return Reference to rOStream.
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
 * @brief Reads an EigenMatrix from the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))".
 * @details Parses the format written by the matching operator<<.
 * @tparam TDataType The matrix's scalar type.
 * @param rIStream Input stream to read from.
 * @param rM Matrix to populate; resized to match the parsed dimensions.
 * @return Reference to rIStream.
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

/**@}*/

} // namespace Kratos

// uBLAS-style free functions (prod, inner_prod, noalias, ...) for the Eigen
// types, kept in a separate header for readability.
#include "includes/eigen_compat_operations.h"
