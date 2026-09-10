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

// See the note in containers/array_1d.h: the Eigen-backed fixed-size types
// rely on the conversions injected into Eigen::MatrixBase by the Kratos
// plugin, which must be active in every translation unit.
#ifndef EIGEN_MATRIXBASE_PLUGIN
#error "The Eigen-backed bounded types require EIGEN_MATRIXBASE_PLUGIN=\"<kratos>/kratos/includes/kratos_eigen_matrixbase_plugin.h\" to be defined globally; configure with KRATOS_LINEAR_ALGEBRA_BACKEND=eigen through the Kratos CMake."
#endif

// System includes
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <type_traits>

// External includes
#include <Eigen/Core>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Internals
{

/// Storage options for the Eigen-backed bounded matrix: row-major to match
/// the uBLAS bounded_matrix memory layout, except where Eigen mandates the
/// orientation (column vectors must be column-major, row vectors row-major).
/// DontAlign keeps the plain T[N1*N2] layout without over-alignment.
constexpr int EigenBoundedMatrixOptions(const std::size_t Size1, const std::size_t Size2)
{
    return ((Size2 == 1 && Size1 != 1) ? Eigen::ColMajor : Eigen::RowMajor) | Eigen::DontAlign;
}

} // namespace Internals

///@name Kratos Classes
///@{

/**
 * @class EigenBoundedMatrix
 * @brief Eigen-backed implementation of the Kratos fixed-size dense matrix (BoundedMatrix).
 * @details Selected by the Eigen linear-algebra backend
 * (KRATOS_USE_EIGEN_BACKEND); under the uBLAS backend BoundedMatrix aliases
 * boost::numeric::ublas::bounded_matrix instead. Mirrors the uBLAS member
 * surface (size1/size2, resize, operator(), the assignment protocol and
 * construction/assignment from uBLAS matrix expressions) so generic Kratos
 * code compiles unchanged. Unlike ublas::bounded_matrix the dimensions are
 * strictly fixed: resizing to anything but (TSize1, TSize2) is an error.
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

    // static bounds of the ublas::bounded_matrix interface (generic code
    // names them, e.g. SVDUtils)
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
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2)))
    EigenBoundedMatrix(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    /// Construction from any uBLAS matrix expression (interoperability with
    /// the dynamic uBLAS Matrix and expressions such as ZeroMatrix or IdentityMatrix).
    template<class TExpression>
    EigenBoundedMatrix(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression) : BaseType()
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size1() != TSize1 || r_expression.size2() != TSize2) << "Wrong size in the construction from a matrix expression [ expression size (" << r_expression.size1() << "," << r_expression.size2() << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
        for (size_type i = 0; i < TSize1; ++i) {
            for (size_type j = 0; j < TSize2; ++j) {
                (*this)(i, j) = r_expression(i, j);
            }
        }
    }

    ///@}
    ///@name Operators
    ///@{

    // Element access
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

    // Assignment
    EigenBoundedMatrix& operator=(const EigenBoundedMatrix& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// The assignment evaluates through a temporary (.eval()), reproducing
    /// the alias-safe uBLAS operator= semantics (e.g. M = trans(M)); the
    /// temporary-free fast path remains noalias(M) = expr, as in uBLAS.
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2)))
    EigenBoundedMatrix& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>().eval());
        return *this;
    }

    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2)))
    EigenBoundedMatrix& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator+=(rExpression.template cast<T>().eval());
        return *this;
    }

    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize1)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(TSize2)))
    EigenBoundedMatrix& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator-=(rExpression.template cast<T>().eval());
        return *this;
    }

    template<class TExpression>
    EigenBoundedMatrix& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size1() != TSize1 || r_expression.size2() != TSize2) << "Wrong size in the assignment from a matrix expression [ expression size (" << r_expression.size1() << "," << r_expression.size2() << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
        for (size_type i = 0; i < TSize1; ++i) {
            for (size_type j = 0; j < TSize2; ++j) {
                (*this)(i, j) = r_expression(i, j);
            }
        }
        return *this;
    }

    template<class TExpression>
    EigenBoundedMatrix& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size1() != TSize1 || r_expression.size2() != TSize2) << "Wrong size in the addition of a matrix expression [ expression size (" << r_expression.size1() << "," << r_expression.size2() << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
        for (size_type i = 0; i < TSize1; ++i) {
            for (size_type j = 0; j < TSize2; ++j) {
                (*this)(i, j) += r_expression(i, j);
            }
        }
        return *this;
    }

    template<class TExpression>
    EigenBoundedMatrix& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size1() != TSize1 || r_expression.size2() != TSize2) << "Wrong size in the subtraction of a matrix expression [ expression size (" << r_expression.size1() << "," << r_expression.size2() << "), bounded size (" << TSize1 << "," << TSize2 << ") ]." << std::endl;
        for (size_type i = 0; i < TSize1; ++i) {
            for (size_type j = 0; j < TSize2; ++j) {
                (*this)(i, j) -= r_expression(i, j);
            }
        }
        return *this;
    }

    // Scalar operators keep the base class implementations (no aliasing there)
    using BaseType::operator*=;
    using BaseType::operator/=;

    /// Implicit conversion to the dynamic uBLAS matrix, mirroring the
    /// implicit matrix_expression conversion ublas::bounded_matrix offers
    /// (e.g. Matrix m = some_bounded_matrix; or DenseVector<Matrix> storage).
    operator boost::numeric::ublas::matrix<T>() const
    {
        boost::numeric::ublas::matrix<T> result(TSize1, TSize2);
        for (size_type i = 0; i < TSize1; ++i) {
            for (size_type j = 0; j < TSize2; ++j) {
                result(i, j) = (*this)(i, j);
            }
        }
        return result;
    }

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
    template<class TExpression>
    EigenBoundedMatrix& assign(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        return (*this = rExpression);
    }
    template<class TExpression>
    EigenBoundedMatrix& plus_assign(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        return (*this += rExpression);
    }
    template<class TExpression>
    EigenBoundedMatrix& minus_assign(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        return (*this -= rExpression);
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

    /// Raw contiguous storage pointer (row-major, matching ublas::bounded_matrix).
    using BaseType::data;

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

    /// static bound of the ublas::bounded_vector interface
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

    /// Construction from any Eigen expression of compatible static size (the
    /// constraint keeps incompatible constructions out of overload sets). The
    /// cast mirrors the uBLAS element-wise conversion semantics for
    /// expressions of a different scalar type.
    template<class TDerived>
    requires ((TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(TSize)) &&
              (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1))
    EigenBoundedVector(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    /// Construction from any uBLAS vector expression.
    template<class TExpression>
    EigenBoundedVector(const boost::numeric::ublas::vector_expression<TExpression>& rExpression) : BaseType()
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != TSize) << "Wrong size in the construction from a vector expression [ expression size = " << r_expression.size() << ", bounded size = " << TSize << " ]." << std::endl;
        for (size_type i = 0; i < TSize; ++i) {
            (*this)[i] = r_expression(i);
        }
    }

    ///@}
    ///@name Operators
    ///@{

    // Element access (both operators)
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

    // Assignment
    EigenBoundedVector& operator=(const EigenBoundedVector& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// The assignment evaluates through a temporary (.eval()), reproducing
    /// the alias-safe uBLAS operator= semantics; the temporary-free fast
    /// path remains noalias(v) = expr, as in uBLAS.
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

    template<class TExpression>
    EigenBoundedVector& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != TSize) << "Wrong size in the assignment from a vector expression [ expression size = " << r_expression.size() << ", bounded size = " << TSize << " ]." << std::endl;
        for (size_type i = 0; i < TSize; ++i) {
            (*this)[i] = r_expression(i);
        }
        return *this;
    }

    template<class TExpression>
    EigenBoundedVector& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != TSize) << "Wrong size in the addition of a vector expression [ expression size = " << r_expression.size() << ", bounded size = " << TSize << " ]." << std::endl;
        for (size_type i = 0; i < TSize; ++i) {
            (*this)[i] += r_expression(i);
        }
        return *this;
    }

    template<class TExpression>
    EigenBoundedVector& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != TSize) << "Wrong size in the subtraction of a vector expression [ expression size = " << r_expression.size() << ", bounded size = " << TSize << " ]." << std::endl;
        for (size_type i = 0; i < TSize; ++i) {
            (*this)[i] -= r_expression(i);
        }
        return *this;
    }

    // Scalar operators keep the base class implementations (no aliasing there)
    using BaseType::operator*=;
    using BaseType::operator/=;

    /// Implicit conversion to the dynamic uBLAS vector, mirroring the
    /// implicit vector_expression conversion ublas::bounded_vector offers.
    operator boost::numeric::ublas::vector<T>() const
    {
        boost::numeric::ublas::vector<T> result(TSize);
        for (size_type i = 0; i < TSize; ++i) {
            result[i] = (*this)[i];
        }
        return result;
    }

    bool operator==(const EigenBoundedVector& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /// Interface-compatibility resize: the size is fixed, resizing to a
    /// different size is an error (unlike ublas::bounded_vector, which can
    /// carry a runtime size below the static maximum).
    void resize(size_type Size, bool Preserve = true)
    {
        KRATOS_ERROR_IF(Size != TSize) << "A bounded vector of static size " << TSize << " cannot be resized to " << Size << " under the Eigen backend." << std::endl;
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    template<class TExpression>
    EigenBoundedVector& assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this = rExpression);
    }
    template<class TExpression>
    EigenBoundedVector& plus_assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this += rExpression);
    }
    template<class TExpression>
    EigenBoundedVector& minus_assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this -= rExpression);
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

    /// Raw contiguous storage pointer.
    using BaseType::data;

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
///@name Input and output
///@{

/// output stream function: reproduces the uBLAS text format [N1,N2]((...),(...))
template<class T, std::size_t TSize1, std::size_t TSize2>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenBoundedMatrix<T, TSize1, TSize2>& rMatrix)
{
    rOStream << '[' << TSize1 << ',' << TSize2 << "](";
    for (std::size_t i = 0; i < TSize1; ++i) {
        if (i > 0) {
            rOStream << ',';
        }
        rOStream << '(';
        for (std::size_t j = 0; j < TSize2; ++j) {
            if (j > 0) {
                rOStream << ',';
            }
            rOStream << rMatrix(i, j);
        }
        rOStream << ')';
    }
    rOStream << ')';
    return rOStream;
}

/// output stream function: reproduces the uBLAS text format [N](v0,v1,...)
template<class T, std::size_t TSize>
inline std::ostream& operator<<(std::ostream& rOStream, const EigenBoundedVector<T, TSize>& rVector)
{
    rOStream << '[' << TSize << "](";
    for (std::size_t i = 0; i < TSize; ++i) {
        if (i > 0) {
            rOStream << ',';
        }
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
