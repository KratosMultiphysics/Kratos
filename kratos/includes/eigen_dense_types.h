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

/// A compile-time vector-shaped Eigen expression (N x 1), the shape uBLAS
/// treats as a vector type.
template<class TDerived>
inline constexpr bool IsColumnShaped = TDerived::ColsAtCompileTime == 1 && TDerived::RowsAtCompileTime != 1;

/// A compile-time vector-shaped Eigen expression of either orientation.
template<class TDerived>
inline constexpr bool IsVectorShaped = TDerived::ColsAtCompileTime == 1 || TDerived::RowsAtCompileTime == 1;

/// Static-shape compatibility of an Eigen expression with a TSize1 x TSize2
/// bounded matrix. Vector shapes of either orientation are compatible with a
/// bounded row/column matrix of the same length (Eigen transposes vectors on
/// assignment), as uBLAS assigns a transposed vector into a 1 x N matrix.
template<class TDerived, std::size_t TSize1, std::size_t TSize2>
inline constexpr bool IsBoundedCompatible =
    ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
     (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2))) ||
    (TSize1 == 1 && TDerived::ColsAtCompileTime == 1 && (TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize2))) ||
    (TSize2 == 1 && TDerived::RowsAtCompileTime == 1 && (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize1)));

/// Storage options for the Eigen-backed bounded matrix: row-major to match
/// the uBLAS bounded_matrix memory layout, except where Eigen mandates the
/// orientation (column vectors must be column-major, row vectors row-major).
/// DontAlign keeps the plain T[N1*N2] layout without over-alignment.
constexpr int EigenBoundedMatrixOptions(const std::size_t Size1, const std::size_t Size2)
{
    return ((Size2 == 1 && Size1 != 1) ? Eigen::ColMajor : Eigen::RowMajor) | Eigen::DontAlign;
}

/// Dynamic, row-major Eigen matrix (the plain object behind the Kratos dense matrices).
template<class T>
using DynamicRowMajorMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/// Dynamic Eigen column vector (the plain object behind the Kratos dense vectors).
template<class T>
using DynamicColumnVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

} // namespace Internals

/**
 * @class EigenMatrix
 * @brief Dense, dynamically sized, row-major Eigen matrix with the uBLAS member surface.
 * @details Row-major storage matches the storage order of the uBLAS dense
 * matrices, so raw-buffer interoperability is preserved. The added members
 * (size1/size2, the preserving resize, clear, data(), the assignment
 * protocol) mirror boost::numeric::ublas::matrix so generic Kratos code
 * compiles unchanged. Eigen only ever sees the base class, so shadowing its
 * members here is safe.
 * @tparam TDataType The scalar type stored in the matrix (e.g. double).
 */
template<class TDataType>
class EigenMatrix : public Internals::DynamicRowMajorMatrix<TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Internals::DynamicRowMajorMatrix<TDataType>;
    using value_type = TDataType;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    using pointer = TDataType*;
    using iterator = TDataType*;
    using const_iterator = const TDataType*;
    using array_type = Internals::StorageView<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: an empty (0 x 0) matrix.
    EigenMatrix() = default;

    EigenMatrix(const EigenMatrix& rOther) = default;
    EigenMatrix(EigenMatrix&& rOther) = default;
    EigenMatrix& operator=(const EigenMatrix& rOther) = default;
    EigenMatrix& operator=(EigenMatrix&& rOther) = default;

    /// Allocates an uninitialized Size1 x Size2 matrix.
    EigenMatrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    /// uBLAS-style (rows, cols, value) fill constructor.
    EigenMatrix(const std::size_t Size1, const std::size_t Size2, const TDataType& rValue) : BaseType(Size1, Size2)
    {
        std::fill(begin(), end(), rValue);
    }

    /// Construction from any matrix-shaped Eigen expression (the scalar cast
    /// makes the cross-precision conversions of ublas work; it is the identity
    /// for matching scalars). Vector-shaped expressions are excluded, as ublas
    /// does not construct a matrix from a vector expression either.
    template<class TDerived>
    requires (!Internals::IsColumnShaped<TDerived>)
    EigenMatrix(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment from any matrix-shaped Eigen expression. The assignment
    /// evaluates through a temporary (.eval()), reproducing the alias-safe
    /// uBLAS operator= semantics (e.g. M = trans(M)); the temporary-free fast
    /// path remains noalias(M) = expr / assign(), as in uBLAS.
    template<class TDerived>
    requires (!Internals::IsColumnShaped<TDerived>)
    EigenMatrix& operator=(const Eigen::MatrixBase<TDerived>& rOther)
    {
        BaseType::operator=(rOther.derived().template cast<TDataType>().eval());
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// uBLAS assignment protocol: assignment without a protective temporary
    /// (resizing the target when needed, a superset of the ublas contract).
    template<class TDerived>
    EigenMatrix& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS assignment protocol: in-place addition without a protective temporary.
    template<class TDerived>
    EigenMatrix& plus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() += rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS assignment protocol: in-place subtraction without a protective temporary.
    template<class TDerived>
    EigenMatrix& minus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() -= rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS-style resize; as in ublas::matrix the default preserves the values.
    void resize(const std::size_t NewSize1, const std::size_t NewSize2, const bool Preserve = true)
    {
        if (Preserve) {
            this->conservativeResize(NewSize1, NewSize2);
        } else {
            BaseType::resize(NewSize1, NewSize2);
        }
    }

    /// uBLAS-style clear: value-initialize every entry (the size is kept).
    void clear()
    {
        std::fill(begin(), end(), TDataType());
    }

    void swap(EigenMatrix& rOther)
    {
        BaseType::swap(rOther);
    }

    friend void swap(EigenMatrix& rFirst, EigenMatrix& rSecond)
    {
        rFirst.swap(rSecond);
    }

    ///@}
    ///@name Access
    ///@{

    /// Number of rows.
    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }

    /// Number of columns.
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }

    /// Contiguous (row-major) storage view, with the uBLAS data() surface and
    /// an implicit conversion to the raw pointer.
    Internals::StorageView<TDataType> data() { return {BaseType::data(), static_cast<std::size_t>(BaseType::size())}; }
    Internals::StorageView<const TDataType> data() const { return {BaseType::data(), static_cast<std::size_t>(BaseType::size())}; }

    iterator begin() { return BaseType::data(); }
    const_iterator begin() const { return BaseType::data(); }
    iterator end() { return BaseType::data() + BaseType::size(); }
    const_iterator end() const { return BaseType::data() + BaseType::size(); }

    ///@}
};

/**
 * @class EigenVector
 * @brief Dense, dynamically sized Eigen column vector with the uBLAS member surface.
 * @details Adds the (size), (size, value) constructors, the preserving resize,
 * the unsigned size(), data() and the assignment protocol of
 * boost::numeric::ublas::vector; operator[] and operator() are already
 * provided by Eigen with compatible semantics.
 * @tparam TDataType The scalar type stored in the vector (e.g. double).
 */
template<class TDataType>
class EigenVector : public Internals::DynamicColumnVector<TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Internals::DynamicColumnVector<TDataType>;
    using value_type = TDataType;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    using pointer = TDataType*;
    using iterator = TDataType*;
    using const_iterator = const TDataType*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using array_type = Internals::StorageView<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: an empty (0-length) vector.
    EigenVector() = default;

    EigenVector(const EigenVector& rOther) = default;
    EigenVector(EigenVector&& rOther) = default;
    EigenVector& operator=(const EigenVector& rOther) = default;
    EigenVector& operator=(EigenVector&& rOther) = default;

    /// Allocates an uninitialized vector of the given size.
    explicit EigenVector(const std::size_t Size) : BaseType(Size) {}

    /// uBLAS-style (size, value) fill constructor.
    EigenVector(const std::size_t Size, const TDataType& rValue) : BaseType(Size)
    {
        std::fill(begin(), end(), rValue);
    }

    /// Construction from any vector-shaped Eigen expression (the scalar cast
    /// makes the cross-precision conversions of ublas work; it is the identity
    /// for matching scalars).
    template<class TDerived>
    requires (Internals::IsVectorShaped<TDerived>)
    EigenVector(const Eigen::MatrixBase<TDerived>& rOther) : BaseType(rOther.derived().template cast<TDataType>()) {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment from any vector-shaped Eigen expression, evaluated through a
    /// temporary (alias-safe, as the ublas operator=); noalias(v) = expr and
    /// assign() remain the temporary-free path.
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

    /// uBLAS assignment protocol: assignment without a protective temporary.
    template<class TDerived>
    EigenVector& assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() = rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS assignment protocol: in-place addition without a protective temporary.
    template<class TDerived>
    EigenVector& plus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() += rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS assignment protocol: in-place subtraction without a protective temporary.
    template<class TDerived>
    EigenVector& minus_assign(const Eigen::MatrixBase<TDerived>& rOther)
    {
        this->noalias() -= rOther.derived().template cast<TDataType>();
        return *this;
    }

    /// uBLAS-style resize; as in ublas::vector the default preserves the values.
    void resize(const std::size_t NewSize, const bool Preserve = true)
    {
        if (Preserve) {
            this->conservativeResize(NewSize);
        } else {
            BaseType::resize(NewSize);
        }
    }

    /// uBLAS-style clear: value-initialize every entry (the size is kept).
    void clear()
    {
        std::fill(begin(), end(), TDataType());
    }

    void swap(EigenVector& rOther)
    {
        BaseType::swap(rOther);
    }

    friend void swap(EigenVector& rFirst, EigenVector& rSecond)
    {
        rFirst.swap(rSecond);
    }

    /// uBLAS-style element insertion (assignment at the given index).
    reference insert_element(const std::size_t Index, const_reference rValue)
    {
        return ((*this)[Index] = rValue);
    }

    /// uBLAS-style element erasure (value-initialization at the given index).
    void erase_element(const std::size_t Index)
    {
        (*this)[Index] = TDataType();
    }

    ///@}
    ///@name Access
    ///@{

    /// uBLAS-style size: unsigned, so the pervasive comparisons/loops against
    /// std::size_t counters stay warning-free (Eigen's own size() is signed).
    std::size_t size() const { return static_cast<std::size_t>(BaseType::size()); }

    /// uBLAS-style emptiness check.
    bool empty() const { return BaseType::size() == 0; }

    /// Contiguous storage view, with the uBLAS data() surface and an implicit
    /// conversion to the raw pointer.
    Internals::StorageView<TDataType> data() { return {BaseType::data(), size()}; }
    Internals::StorageView<const TDataType> data() const { return {BaseType::data(), size()}; }

    iterator begin() { return BaseType::data(); }
    const_iterator begin() const { return BaseType::data(); }
    const_iterator cbegin() const { return BaseType::data(); }
    iterator end() { return BaseType::data() + BaseType::size(); }
    const_iterator end() const { return BaseType::data() + BaseType::size(); }
    const_iterator cend() const { return BaseType::data() + BaseType::size(); }

    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    ///@}
};

/**
 * @class EigenBoundedMatrix
 * @brief Eigen-backed implementation of the Kratos fixed-size dense matrix (BoundedMatrix).
 * @details Mirrors the uBLAS bounded_matrix member surface (size1/size2,
 * resize, operator(), the assignment protocol, iterators) so generic Kratos
 * code compiles unchanged. Unlike ublas::bounded_matrix the dimensions are
 * strictly fixed: resizing to anything but (TSize1, TSize2) is an error (use
 * a dynamic Matrix where a variable size is needed).
 */
template<class T, std::size_t TSize1, std::size_t TSize2>
class EigenBoundedMatrix : public Eigen::Matrix<T, static_cast<int>(TSize1), static_cast<int>(TSize2), Internals::EigenBoundedMatrixOptions(TSize1, TSize2)>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigenBoundedMatrix);

    using BaseType = Eigen::Matrix<T, static_cast<int>(TSize1), static_cast<int>(TSize2), Internals::EigenBoundedMatrixOptions(TSize1, TSize2)>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using const_reference = const T&;
    using reference = T&;
    using pointer = T*;
    using iterator = T*;
    using const_iterator = const T*;
    using array_type = Internals::StorageView<T>;

    /// Static bounds of the ublas::bounded_matrix interface (named by generic code, e.g. SVDUtils).
    static constexpr size_type max_size1 = TSize1;
    static constexpr size_type max_size2 = TSize2;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: full (TSize1 x TSize2) size, values NOT initialized.
    EigenBoundedMatrix() : BaseType() {}

    /// The size arguments are only there for interface compatibility.
    EigenBoundedMatrix(size_type Size1, size_type Size2) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size1 != TSize1 || Size2 != TSize2) << "Wrong size in the construction of a bounded matrix [ requested (" << Size1 << "," << Size2 << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
    }

    EigenBoundedMatrix(const EigenBoundedMatrix& rOther) : BaseType(rOther) {}

    /// Construction from any Eigen expression of compatible static size (the
    /// constraint keeps incompatible constructions out of overload sets). The
    /// cast mirrors the uBLAS element-wise conversion semantics for
    /// expressions of a different scalar type.
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    ///@}
    ///@name Operators
    ///@{

    const_reference operator()(size_type I, size_type J) const
    {
        KRATOS_DEBUG_ERROR_IF(I >= TSize1 || J >= TSize2) << "Index greater than the size of the bounded matrix - indices are (" << I << "," << J << ")" << std::endl;
        return BaseType::coeff(I, J);
    }
    reference operator()(size_type I, size_type J)
    {
        KRATOS_DEBUG_ERROR_IF(I >= TSize1 || J >= TSize2) << "Index greater than the size of the bounded matrix - indices are (" << I << "," << J << ")" << std::endl;
        return BaseType::coeffRef(I, J);
    }

    EigenBoundedMatrix& operator=(const EigenBoundedMatrix& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// The assignment evaluates through a temporary (.eval()), reproducing
    /// the alias-safe uBLAS operator= semantics (e.g. M = trans(M)); the
    /// temporary-free fast path remains noalias(M) = expr, as in uBLAS.
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>().eval());
        return *this;
    }

    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, TSize1, TSize2>)
    EigenBoundedMatrix& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator+=(rExpression.template cast<T>().eval());
        return *this;
    }

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

    bool operator==(const EigenBoundedMatrix& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /// Interface-compatibility resize: the dimensions are fixed, resizing to
    /// a different size is an error (unlike ublas::bounded_matrix, which can
    /// carry runtime sizes below the static maximum).
    void resize(size_type Size1, size_type Size2, bool Preserve = true)
    {
        KRATOS_ERROR_IF(Size1 != TSize1 || Size2 != TSize2) << "A bounded matrix of static size (" << TSize1 << "," << TSize2 << ") cannot be resized to (" << Size1 << "," << Size2 << ") under the Eigen backend." << std::endl;
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    template<class TDerived>
    EigenBoundedMatrix& assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() = rExpression.template cast<T>();
        return *this;
    }
    template<class TDerived>
    EigenBoundedMatrix& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() += rExpression.template cast<T>();
        return *this;
    }
    template<class TDerived>
    EigenBoundedMatrix& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() -= rExpression.template cast<T>();
        return *this;
    }

    void swap(EigenBoundedMatrix& rOther)
    {
        if (this != &rOther) {
            std::swap_ranges(begin(), end(), rOther.begin());
        }
    }
    friend void swap(EigenBoundedMatrix& rFirst, EigenBoundedMatrix& rSecond)
    {
        rFirst.swap(rSecond);
    }

    void clear()
    {
        std::fill(begin(), end(), value_type());
    }

    ///@}
    ///@name Access
    ///@{

    constexpr size_type size1() const { return TSize1; }
    constexpr size_type size2() const { return TSize2; }

    /// Contiguous (row-major) storage view, with the uBLAS data() surface and
    /// an implicit conversion to the raw pointer.
    Internals::StorageView<T> data() { return {BaseType::data(), TSize1 * TSize2}; }
    Internals::StorageView<const T> data() const { return {BaseType::data(), TSize1 * TSize2}; }

    iterator begin() { return BaseType::data(); }
    const_iterator begin() const { return BaseType::data(); }
    iterator end() { return BaseType::data() + TSize1 * TSize2; }
    const_iterator end() const { return BaseType::data() + TSize1 * TSize2; }

    ///@}
}; // Class EigenBoundedMatrix

/**
 * @class EigenBoundedVector
 * @brief Eigen-backed implementation of the Kratos fixed-size dense vector alias (BoundedVector).
 * @details Same surface as the Eigen array_1d but a distinct type, since
 * generic code (e.g. the serializer) overloads on array_1d and BoundedVector
 * simultaneously. Unlike ublas::bounded_vector the size is strictly fixed:
 * resizing to anything but TSize is an error.
 */
template<class T, std::size_t TSize>
class EigenBoundedVector : public Eigen::Matrix<T, static_cast<int>(TSize), 1, Eigen::ColMajor | Eigen::DontAlign>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigenBoundedVector);

    using BaseType = Eigen::Matrix<T, static_cast<int>(TSize), 1, Eigen::ColMajor | Eigen::DontAlign>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using const_reference = const T&;
    using reference = T&;
    using pointer = T*;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using array_type = Internals::StorageView<T>;

    /// Static bound of the ublas::bounded_vector interface.
    static constexpr size_type max_size = TSize;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: full size, values NOT initialized.
    EigenBoundedVector() : BaseType() {}

    /// The size argument is only there for interface compatibility.
    explicit EigenBoundedVector(size_type Size) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size != TSize) << "Wrong size in the construction of a bounded vector [ requested " << Size << ", bounded size " << TSize << " ]." << std::endl;
    }

    EigenBoundedVector(size_type Size, const value_type& rValue) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size > TSize) << "Given size is greater than the size of the bounded vector!" << std::endl;
        std::fill(begin(), begin() + Size, rValue);
    }

    EigenBoundedVector(const EigenBoundedVector& rOther) : BaseType(rOther) {}

    /// Construction from any Eigen expression of compatible static size.
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    ///@}
    ///@name Operators
    ///@{

    const_reference operator()(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    reference operator()(size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    const_reference operator[](size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    reference operator[](size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= TSize) << "Index greater than the size of the bounded vector - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    EigenBoundedVector& operator=(const EigenBoundedVector& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// Alias-safe assignment (through a temporary), as the ublas operator=.
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>().eval());
        return *this;
    }

    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator+=(rExpression.template cast<T>().eval());
        return *this;
    }

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

    bool operator==(const EigenBoundedVector& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /// Interface-compatibility resize: the size is fixed, resizing to a
    /// different size is an error (unlike ublas::bounded_vector).
    void resize(size_type Size, bool Preserve = true)
    {
        KRATOS_ERROR_IF(Size != TSize) << "A bounded vector of static size " << TSize << " cannot be resized to " << Size << " under the Eigen backend." << std::endl;
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    template<class TDerived>
    EigenBoundedVector& assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() = rExpression.template cast<T>();
        return *this;
    }
    template<class TDerived>
    EigenBoundedVector& plus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() += rExpression.template cast<T>();
        return *this;
    }
    template<class TDerived>
    EigenBoundedVector& minus_assign(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        this->noalias() -= rExpression.template cast<T>();
        return *this;
    }

    void swap(EigenBoundedVector& rOther)
    {
        if (this != &rOther) {
            std::swap_ranges(begin(), end(), rOther.begin());
        }
    }
    friend void swap(EigenBoundedVector& rFirst, EigenBoundedVector& rSecond)
    {
        rFirst.swap(rSecond);
    }

    void clear()
    {
        std::fill(begin(), end(), value_type());
    }

    ///@}
    ///@name Access
    ///@{

    constexpr size_type size() const
    {
        return TSize;
    }

    /// Contiguous storage view, with the uBLAS data() surface and an implicit
    /// conversion to the raw pointer.
    Internals::StorageView<T> data() { return {BaseType::data(), TSize}; }
    Internals::StorageView<const T> data() const { return {BaseType::data(), TSize}; }

    iterator begin() { return BaseType::data(); }
    const_iterator begin() const { return BaseType::data(); }
    const_iterator cbegin() const { return BaseType::data(); }
    iterator end() { return BaseType::data() + TSize; }
    const_iterator end() const { return BaseType::data() + TSize; }
    const_iterator cend() const { return BaseType::data() + TSize; }

    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
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
/// Nullary functor of the unit vector e_i.
template<class T>
struct UnitVectorOp
{
    Eigen::Index mIndex;
    T operator()(const Eigen::Index Index) const { return Index == mIndex ? T(1) : T(0); }
};
} // namespace Internals

/// Lazy zero matrix (uBLAS zero_matrix): zero_matrix(n) is n x n, zero_matrix(n, m) is n x m.
template<class T = double>
class zero_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>;
    using value_type = T;
    using size_type = std::size_t;

    zero_matrix(const zero_matrix& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    zero_matrix& operator=(const zero_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    zero_matrix() : BaseType(0, 0, Eigen::internal::scalar_constant_op<T>(T(0))) {}
    explicit zero_matrix(const std::size_t Size) : BaseType(Size, Size, Eigen::internal::scalar_constant_op<T>(T(0))) {}
    zero_matrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/// Lazy identity matrix (uBLAS identity_matrix): identity_matrix(n) is n x n, identity_matrix(n, m) is n x m.
template<class T = double>
class identity_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<T>, Internals::DynamicRowMajorMatrix<T>>;
    using value_type = T;
    using size_type = std::size_t;

    identity_matrix(const identity_matrix& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    identity_matrix& operator=(const identity_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    identity_matrix() : BaseType(0, 0) {}
    explicit identity_matrix(const std::size_t Size) : BaseType(Size, Size) {}
    identity_matrix(const std::size_t Size1, const std::size_t Size2) : BaseType(Size1, Size2) {}

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/// Lazy constant matrix (uBLAS scalar_matrix): every entry equals the given value (1 by default, as in uBLAS).
template<class T = double>
class scalar_matrix : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicRowMajorMatrix<T>>;
    using value_type = T;
    using size_type = std::size_t;

    scalar_matrix(const scalar_matrix& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    scalar_matrix& operator=(const scalar_matrix& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    scalar_matrix() : BaseType(0, 0, Eigen::internal::scalar_constant_op<T>(T(1))) {}
    scalar_matrix(const std::size_t Size1, const std::size_t Size2, const T& rValue = T(1)) : BaseType(Size1, Size2, Eigen::internal::scalar_constant_op<T>(rValue)) {}

    std::size_t size1() const { return static_cast<std::size_t>(this->rows()); }
    std::size_t size2() const { return static_cast<std::size_t>(this->cols()); }
};

/// Lazy zero vector (uBLAS zero_vector).
template<class T = double>
class zero_vector : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>;
    using value_type = T;
    using size_type = std::size_t;

    zero_vector(const zero_vector& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    zero_vector& operator=(const zero_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    zero_vector() : BaseType(0, 1, Eigen::internal::scalar_constant_op<T>(T(0))) {}
    explicit zero_vector(const std::size_t Size) : BaseType(Size, 1, Eigen::internal::scalar_constant_op<T>(T(0))) {}

    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }
};

/// Lazy constant vector (uBLAS scalar_vector): every entry equals the given value (1 by default, as in uBLAS).
template<class T = double>
class scalar_vector : public Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<T>, Internals::DynamicColumnVector<T>>;
    using value_type = T;
    using size_type = std::size_t;

    scalar_vector(const scalar_vector& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    scalar_vector& operator=(const scalar_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    scalar_vector() : BaseType(0, 1, Eigen::internal::scalar_constant_op<T>(T(1))) {}
    explicit scalar_vector(const std::size_t Size, const T& rValue = T(1)) : BaseType(Size, 1, Eigen::internal::scalar_constant_op<T>(rValue)) {}

    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }
};

/// Lazy unit vector e_i (uBLAS unit_vector).
template<class T = double>
class unit_vector : public Eigen::CwiseNullaryOp<Internals::UnitVectorOp<T>, Internals::DynamicColumnVector<T>>
{
public:
    using BaseType = Eigen::CwiseNullaryOp<Internals::UnitVectorOp<T>, Internals::DynamicColumnVector<T>>;
    using value_type = T;
    using size_type = std::size_t;

    unit_vector(const unit_vector& rOther) = default;
    /// Assignable like the uBLAS factories: the Eigen nullary base has const
    /// extents and functor, so the value object is rebuilt in place.
    unit_vector& operator=(const unit_vector& rOther)
    {
        if (this != &rOther) {
            std::destroy_at(this);
            std::construct_at(this, rOther);
        }
        return *this;
    }

    unit_vector() : BaseType(0, 1, Internals::UnitVectorOp<T>{0}) {}
    explicit unit_vector(const std::size_t Size, const std::size_t Index = 0) : BaseType(Size, 1, Internals::UnitVectorOp<T>{static_cast<Eigen::Index>(Index)}) {}

    std::size_t size() const { return static_cast<std::size_t>(this->rows()); }
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

/// Writes an EigenVector in the uBLAS text format "[n](v0,v1,...)".
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

/// Reads an EigenVector from the uBLAS text format "[n](v0,v1,...)"; the
/// failbit is set on any mismatch, as boost's own operator>> does.
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

/// Writes an EigenMatrix in the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))".
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

/// Reads an EigenMatrix from the uBLAS text format "[n1,n2]((a00,a01),(a10,a11))".
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

/// Writes an EigenBoundedMatrix in the uBLAS text format "[N1,N2]((...),(...))".
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

/// Writes an EigenBoundedVector in the uBLAS text format "[N](v0,v1,...)".
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
