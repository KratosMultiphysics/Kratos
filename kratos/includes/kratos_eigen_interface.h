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
#include <type_traits>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "includes/exception.h"

namespace Kratos {

/// Index type used for the Eigen sparse matrices. Eigen requires a signed
/// StorageIndex; std::ptrdiff_t matches the width of the std::size_t indices
/// used by uBLAS and the DOF equation ids, so no narrowing can occur. It can
/// be overridden at configure time with -DKRATOS_EIGEN_INDEX_TYPE=<type>.
#ifdef KRATOS_EIGEN_INDEX_TYPE
using KratosEigenIndexType = KRATOS_EIGEN_INDEX_TYPE;
#else
using KratosEigenIndexType = std::ptrdiff_t;
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
 */
template<class T>
class EigenArrayProxy
{
public:
    using value_type = std::remove_const_t<T>;
    using iterator = T*;
    using const_iterator = const T*;

    EigenArrayProxy(T* pData, const std::size_t Size) : mpData(pData), mSize(Size) {}

    T* begin() { return mpData; }
    T* end() { return mpData + mSize; }
    const T* begin() const { return mpData; }
    const T* end() const { return mpData + mSize; }

    T& operator[](const std::size_t Index) { return mpData[Index]; }
    const T& operator[](const std::size_t Index) const { return mpData[Index]; }

    std::size_t size() const { return mSize; }

private:
    T* mpData;
    std::size_t mSize;
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
 */
template<class TDataType>
class EigenMatrix : public Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
public:
    using BaseType = Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using value_type = TDataType;
    using size_type = std::size_t;

    EigenMatrix() = default;

    EigenMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /// Construction from any Eigen expression
    template<class TDerived>
    EigenMatrix(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /// Assignment from any Eigen expression
    template<class TDerived>
    EigenMatrix& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /// uBLAS-style resize; note that, as in ublas::matrix, the default preserves the values.
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize1, NewSize2);
        else BaseType::resize(NewSize1, NewSize2);
    }
};

/**
 * @class EigenVector
 * @brief Dense, dynamically sized Eigen column vector with the uBLAS member surface.
 * @details Adds the (size), (size, value) constructors and the 2-argument
 * resize of boost::numeric::ublas::vector; operator[], operator() and size()
 * are already provided by Eigen with compatible semantics.
 */
template<class TDataType>
class EigenVector : public Eigen::Matrix<TDataType, Eigen::Dynamic, 1>
{
public:
    using BaseType = Eigen::Matrix<TDataType, Eigen::Dynamic, 1>;
    using value_type = TDataType;
    using size_type = std::size_t;

    EigenVector() = default;

    explicit EigenVector(const std::size_t Size) : BaseType(Size) {}

    EigenVector(const std::size_t Size, const TDataType Value)
        : BaseType(BaseType::Constant(Size, Value)) {}

    /// Construction from any Eigen expression
    template<class TDerived>
    EigenVector(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /// Assignment from any Eigen expression
    template<class TDerived>
    EigenVector& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// uBLAS-style resize; note that, as in ublas::vector, the default preserves the values.
    void resize(const std::size_t NewSize, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize);
        else BaseType::resize(NewSize);
    }

    /// uBLAS-style clear: zero all entries (the size is kept).
    void clear()
    {
        this->setZero();
    }
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
 * @warning As for the uBLAS type, element insertion through operator() on a
 * missing entry is an O(nnz-in-row) slow path and must not be used for assembly.
 */
template<class TDataType, class TIndexType = KratosEigenIndexType>
class EigenCompressedMatrix : public Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>
{
public:
    using BaseType = Eigen::SparseMatrix<TDataType, Eigen::RowMajor, TIndexType>;
    using value_type = TDataType;
    using size_type = std::size_t;
    // uBLAS-style storage-array typedefs (generic code uses e.g.
    // typename TMatrix::index_array_type::value_type to name the index type)
    using index_array_type = Internals::EigenArrayProxy<TIndexType>;
    using value_array_type = Internals::EigenArrayProxy<TDataType>;

    EigenCompressedMatrix() = default;

    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /// uBLAS-style (rows, cols, nnz) constructor: allocates the compressed
    /// storage up front so the CSR arrays can be written directly.
    EigenCompressedMatrix(const std::size_t Size1, const std::size_t Size2, const std::size_t NNZ)
        : BaseType(Size1, Size2)
    {
        this->resizeNonZeros(NNZ);
    }

    /// Construction from any Eigen sparse expression
    template<class TDerived>
    EigenCompressedMatrix(const Eigen::SparseMatrixBase<TDerived>& rOther) : BaseType(rOther) {}

    /// Assignment from any Eigen sparse expression
    template<class TDerived>
    EigenCompressedMatrix& operator=(const Eigen::SparseMatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /// Number of stored entries. As in ublas::compressed_matrix this is the
    /// *filled* count (in Eigen's compressed mode it is derived from the row
    /// pointers), which is 0 until the CSR arrays are written; the allocated
    /// capacity is exposed through the size of the value/index proxies.
    std::size_t nnz() const { return static_cast<std::size_t>(this->nonZeros()); }

    // The storage proxies span the *allocated* capacity (as the uBLAS storage
    // arrays do), so they can be written before the row pointers declare the
    // filled size.
    auto value_data() { return Internals::EigenArrayProxy<TDataType>(this->valuePtr(), this->data().size()); }
    auto index1_data() { return Internals::EigenArrayProxy<TIndexType>(this->outerIndexPtr(), size1() + 1); }
    auto index2_data() { return Internals::EigenArrayProxy<TIndexType>(this->innerIndexPtr(), this->data().size()); }
    auto value_data() const { return Internals::EigenArrayProxy<const TDataType>(this->valuePtr(), this->data().size()); }
    auto index1_data() const { return Internals::EigenArrayProxy<const TIndexType>(this->outerIndexPtr(), size1() + 1); }
    auto index2_data() const { return Internals::EigenArrayProxy<const TIndexType>(this->innerIndexPtr(), this->data().size()); }

    /// uBLAS-style finalization after writing the CSR arrays by hand. The
    /// storage was already sized by the (rows, cols, nnz) constructor and the
    /// filled count is implied by the written row pointers, so this only
    /// validates consistency.
    void set_filled(const std::size_t FilledSize1, const std::size_t FilledNNZ)
    {
        KRATOS_DEBUG_ERROR_IF(FilledSize1 != size1() + 1 || FilledNNZ != nnz() || FilledNNZ > this->data().size())
            << "set_filled(" << FilledSize1 << ", " << FilledNNZ
            << ") is inconsistent with the written compressed storage ("
            << size1() + 1 << ", " << nnz() << " of " << this->data().size()
            << " allocated)." << std::endl;
    }

    /// uBLAS-style resize; as in ublas::compressed_matrix the default preserves,
    /// while resize(m, n, false) discards both values and structure.
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize1, NewSize2);
        else BaseType::resize(NewSize1, NewSize2);
    }

    /// uBLAS-style clear: removes all stored entries (the dimensions are kept).
    void clear()
    {
        this->setZero();
    }

    /// uBLAS-style element insertion. As for operator() on a missing entry
    /// this is an O(nnz-in-row) slow path meant for tests and small setup
    /// code, not for assembly (write the CSR arrays directly instead).
    void push_back(const std::size_t I, const std::size_t J, const TDataType Value)
    {
        this->coeffRef(I, J) = Value;
    }

    /// uBLAS-style element access. Inserting a new entry through the non-const
    /// version is the same slow path as for ublas::compressed_matrix.
    TDataType& operator()(const std::size_t I, const std::size_t J) { return this->coeffRef(I, J); }
    TDataType operator()(const std::size_t I, const std::size_t J) const { return this->coeff(I, J); }
};

} // namespace Kratos

// uBLAS-style free functions (prod, inner_prod, noalias, ...) for the Eigen
// types, kept in a separate header for readability.
#include "includes/eigen_compat_operations.h"
