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

/// Index type used for the Eigen sparse matrices. Eigen requires a signed
/// StorageIndex. The default is a 32-bit int: sparse kernels are memory-bound
/// (SpMV traffic drops from 16 to 12 bytes per nonzero versus 64-bit indices,
/// a structural advantage the std::size_t-indexed uBLAS backend cannot match),
/// and 2^31 - 1 nonzeros per matrix (~16 GB of values alone) comfortably
/// covers shared-memory problems. For larger systems configure with
/// KRATOS_EIGEN_64BIT_INDICES=ON (or override the type directly with
/// -DKRATOS_EIGEN_INDEX_TYPE=<type>).
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

    // Minimal uBLAS container-trait surface, mirroring the EigenVector one:
    // the generic boost::numeric::ublas::project() overloads instantiate
    // matrix_range<M>/matrix_slice<M> while they are *considered* during
    // overload resolution, and that class instantiation needs these member
    // typedefs to be well-formed. The proxies are never actually used on this
    // type - the Kratos Eigen overloads win the resolution - so the closure
    // aliases are plain references and the (never dereferenced through the
    // uBLAS iterator protocol) iterator aliases are plain pointers.
    using difference_type = std::ptrdiff_t;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    using closure_type = EigenMatrix<TDataType>&;
    using const_closure_type = const EigenMatrix<TDataType>&;
    using storage_category = boost::numeric::ublas::dense_tag;
    using orientation_category = boost::numeric::ublas::row_major_tag;
    using iterator1 = TDataType*;
    using const_iterator1 = const TDataType*;
    using iterator2 = TDataType*;
    using const_iterator2 = const TDataType*;

    EigenMatrix() = default;

    EigenMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /// uBLAS-style (rows, cols, value) fill constructor
    EigenMatrix(const std::size_t Size1, const std::size_t Size2, const TDataType Value)
        : BaseType(BaseType::Constant(Size1, Size2, Value)) {}

    /// Construction from any Eigen expression. The scalar cast makes the
    /// cross-precision conversions of ublas (double <-> float grids in the
    /// p-multigrid hierarchy) work; it is the identity for matching scalars.
    template<class TDerived>
    EigenMatrix(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    /// Assignment from any Eigen expression (scalar-casting, see the constructor)
    template<class TDerived>
    EigenMatrix& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>());
        return *this;
    }

    /// Construction from a dense uBLAS matrix expression (ZeroMatrix,
    /// IdentityMatrix, prod(...) results, ...), mirroring ublas::matrix's own
    /// converting constructor so uBLAS idioms keep compiling on this type.
    template<class TExpression>
    EigenMatrix(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
        : BaseType(rExpression().size1(), rExpression().size2())
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size1(); ++i)
            for (std::size_t j = 0; j < r_e.size2(); ++j)
                (*this)(i, j) = r_e(i, j);
    }

    /// Assignment from a dense uBLAS matrix expression (resizing, as ublas does)
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

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /// uBLAS-style resize; note that, as in ublas::matrix, the default preserves the values.
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) this->conservativeResize(NewSize1, NewSize2);
        else BaseType::resize(NewSize1, NewSize2);
    }

    /// uBLAS-style clear: zero all entries (the size is kept).
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
 */
template<class TDataType>
class EigenVector : public Eigen::Matrix<TDataType, Eigen::Dynamic, 1>
{
public:
    using BaseType = Eigen::Matrix<TDataType, Eigen::Dynamic, 1>;
    using value_type = TDataType;
    using size_type = std::size_t;

    // Minimal uBLAS container-trait surface. The generic
    // boost::numeric::ublas::project()/subrange() overloads instantiate
    // vector_range<V>/vector_slice<V> while they are *considered* during
    // overload resolution (their second parameter names a nested type of the
    // proxy), and that class instantiation needs these member typedefs to be
    // well-formed. The proxies are never actually used on this type - the
    // Kratos Eigen overloads win the resolution - so reference-typed closure
    // aliases are sufficient.
    using difference_type = std::ptrdiff_t;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    // (iterator/const_iterator are the inherited Eigen STL iterator typedefs.)
    using closure_type = EigenVector<TDataType>&;
    using const_closure_type = const EigenVector<TDataType>&;
    using storage_category = boost::numeric::ublas::dense_tag;

    EigenVector() = default;

    explicit EigenVector(const std::size_t Size) : BaseType(Size) {}

    EigenVector(const std::size_t Size, const TDataType Value)
        : BaseType(BaseType::Constant(Size, Value)) {}

    /// Construction from any Eigen expression. The scalar cast makes the
    /// cross-precision conversions of ublas (double <-> float grids in the
    /// p-multigrid hierarchy) work; it is the identity for matching scalars.
    template<class TDerived>
    EigenVector(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    /// Assignment from any Eigen expression (scalar-casting, see the constructor)
    template<class TDerived>
    EigenVector& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>());
        return *this;
    }

    /// Construction from a dense uBLAS vector expression (ZeroVector,
    /// prod(...) results, ...), mirroring ublas::vector's own converting
    /// constructor so uBLAS idioms keep compiling on this type.
    template<class TExpression>
    EigenVector(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
        : BaseType(rExpression().size())
    {
        const auto& r_e = rExpression();
        for (std::size_t i = 0; i < r_e.size(); ++i) (*this)[i] = r_e(i);
    }

    /// Assignment from a dense uBLAS vector expression (resizing, as ublas does)
    template<class TExpression>
    EigenVector& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        BaseType::resize(r_e.size());
        for (std::size_t i = 0; i < r_e.size(); ++i) (*this)[i] = r_e(i);
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

    /// uBLAS-style size: unsigned, so the pervasive comparisons/loops against
    /// std::size_t counters stay warning-free (Eigen's own size() is signed).
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }

    /// uBLAS-style emptiness check
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

    /// Construction from a dense Eigen expression: the nonzero entries are
    /// gathered into compressed storage (the dense counterpart of the sparse
    /// constructor above), mirroring ublas::compressed_matrix's converting
    /// constructor for the Eigen-backed dense types.
    template<class TDerived>
    explicit EigenCompressedMatrix(const Eigen::MatrixBase<TDerived>& rExpression)
        : BaseType(rExpression.derived().sparseView())
    {
        this->makeCompressed();
    }

    /// Construction from a dense uBLAS matrix expression, mirroring
    /// ublas::compressed_matrix's own constructor: the nonzero entries are
    /// gathered into compressed storage. Insertion happens in row-major order
    /// (this matrix's storage order), so it is a single linear fill. Not
    /// explicit: EigenMatrix's equivalent constructor isn't either, and the
    /// ublas idiom `SparseMatrixType m = ZeroMatrix(n, n);` relies on the
    /// implicit conversion (as does its uBLAS-backend counterpart).
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
        KRATOS_DEBUG_ERROR_IF(FilledSize1 != size1() + 1 || FilledNNZ != nnz() || FilledNNZ > static_cast<std::size_t>(this->data().size()))
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

///@name uBLAS-format stream operators
///@{
// The dynamic Eigen-backed types print and parse in the boost::numeric::ublas
// text format ("[n](v0,v1,...)" and "[n1,n2]((a00,a01),(a10,a11))"): the
// format is part of the established Kratos IO surface (mdpa Matrix/Vector
// values, json checks, restart files, printed expectations in tests), so it
// must not change with the backend. These exact-type overloads win over
// Eigen's generic MatrixBase operator<<.

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

template<class TDataType>
inline std::istream& operator>>(std::istream& rIStream, EigenVector<TDataType>& rV)
{
    // Parses the ublas vector format written above. On any mismatch the
    // failbit is set, as boost's own operator>> does.
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

template<class TDataType>
inline std::istream& operator>>(std::istream& rIStream, EigenMatrix<TDataType>& rM)
{
    // Parses the ublas matrix format written above.
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

///@}

} // namespace Kratos

// uBLAS-style free functions (prod, inner_prod, noalias, ...) for the Eigen
// types, kept in a separate header for readability.
#include "includes/eigen_compat_operations.h"
