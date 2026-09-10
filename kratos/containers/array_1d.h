//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#pragma once

// The fixed-size dense vector of Kratos. The implementation follows the
// linear-algebra backend selected at configure time
// (KRATOS_LINEAR_ALGEBRA_BACKEND): an Eigen fixed-size vector under the Eigen
// backend, the historical boost::numeric::ublas vector_expression otherwise.
// Both implementations honour the same contract (see
// tests/cpp_tests/containers/test_dense_backend_contract.cpp): plain T[N]
// layout starting at the object address, the uBLAS-style member surface
// (operator()/operator[], size, resize, iterators, ...) and the uBLAS text
// format for streaming.

#ifdef KRATOS_USE_EIGEN_BACKEND

// The Kratos CMake defines EIGEN_MATRIXBASE_PLUGIN globally when the Eigen
// backend is selected; the injected conversions (Eigen expression -> uBLAS
// dense containers) are part of the array_1d contract, so refuse to build
// without them (e.g. a hand-rolled build defining only the backend macro).
#ifndef EIGEN_MATRIXBASE_PLUGIN
#error "KRATOS_USE_EIGEN_BACKEND requires EIGEN_MATRIXBASE_PLUGIN=\"<kratos>/kratos/includes/kratos_eigen_matrixbase_plugin.h\" to be defined globally; configure with KRATOS_LINEAR_ALGEBRA_BACKEND=eigen through the Kratos CMake."
#endif

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <ostream>
#include <type_traits>

// External includes
#include <Eigen/Core>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/vector.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h" // also pulls in the Eigen compatibility operations under this backend

namespace Kratos
{

namespace Internals
{
/// An Eigen expression whose static shape is compatible with an N x 1 column
/// vector (the array_1d layout).
template<class TDerived, std::size_t N>
inline constexpr bool IsColumnShapedEigenExpression =
    (TDerived::RowsAtCompileTime == Eigen::Dynamic || TDerived::RowsAtCompileTime == static_cast<int>(N)) &&
    (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == 1);

/// A statically row-shaped (1 x N) Eigen expression: the uBLAS world treats a
/// matrix row as a plain vector, so array_1d accepts it transposed.
template<class TDerived, std::size_t N>
inline constexpr bool IsRowShapedEigenExpression =
    !IsColumnShapedEigenExpression<TDerived, N> &&
    TDerived::RowsAtCompileTime == 1 &&
    (TDerived::ColsAtCompileTime == Eigen::Dynamic || TDerived::ColsAtCompileTime == static_cast<int>(N));
} // namespace Internals

///@name Kratos Classes
///@{

/**
 * @class array_1d
 * @brief Eigen-backed implementation of the Kratos fixed-size dense vector.
 * @details Selected by the Eigen linear-algebra backend
 * (KRATOS_USE_EIGEN_BACKEND); the uBLAS backend uses the historical
 * boost::numeric::ublas vector_expression implementation instead (below).
 * Both implementations honour the same contract:
 * - Plain T[N] memory layout starting at the object address (no vtable, no
 *   padding, no over-alignment). Variable component access, DataTypeTraits
 *   contiguity and the MPI buffers rely on it, hence Eigen::DontAlign and the
 *   static_asserts at the end of this block.
 * - The uBLAS-style member surface: operator() and operator[], size(),
 *   near-no-op resize(), clear(), raw-pointer iterators, swap, the
 *   assign/plus_assign/minus_assign protocol and construction/assignment from
 *   uBLAS vector expressions (the dynamic Vector/Matrix types remain uBLAS
 *   under both backends).
 * - The uBLAS text format for streaming: [N](v0,v1,...).
 * Deriving from Eigen::Matrix makes every Eigen expression (a+b, prod, ...)
 * available on the type; the mixed uBLAS/Eigen idiom layer lives in
 * includes/eigen_ublas_compat_operations.h.
 */
template<class T, std::size_t N>
class array_1d : public Eigen::Matrix<T, static_cast<int>(N), 1, Eigen::ColMajor | Eigen::DontAlign>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of array_1d
    KRATOS_CLASS_POINTER_DEFINITION(array_1d);

    using BaseType = Eigen::Matrix<T, static_cast<int>(N), 1, Eigen::ColMajor | Eigen::DontAlign>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using const_reference = const T&;
    using reference = T&;
    using array_type = std::array<T, N>;
    using pointer = T*;
    using self_type = array_1d<T, N>;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor: as in the uBLAS implementation the values are NOT initialized.
    array_1d() : BaseType() {}

    /// The size argument is only there for interface compatibility: the size is fixed to N.
    explicit array_1d(size_type Size) : BaseType() {}

    /// Fills the first Size entries with the given value.
    explicit array_1d(size_type Size, const value_type& rValue) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(Size > N) << "Given size is greater than the size of the array!" << std::endl;
        std::fill(begin(), begin() + Size, rValue);
    }

    explicit array_1d(const std::initializer_list<value_type>& rInitList) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(rInitList.size() > N) << "Size of list greater than the size of the array!" << std::endl;
        std::copy(rInitList.begin(), rInitList.end(), begin());
    }

    array_1d(size_type Size, const array_type& rData) : BaseType()
    {
        std::copy(rData.begin(), rData.end(), begin());
    }

    array_1d(const array_1d& rOther) : BaseType(rOther) {}

    /// Construction from any Eigen expression of compatible static size (the
    /// constraint keeps incompatible constructions out of overload sets and
    /// turns size bugs into clear no-viable-overload errors). The cast
    /// mirrors the uBLAS element-wise conversion semantics for expressions of
    /// a different scalar type (e.g. array_1d<int, 3> from an
    /// array_1d<std::size_t, 3>).
    template<class TDerived>
    requires (Internals::IsColumnShapedEigenExpression<TDerived, N>)
    array_1d(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

    /// Row-shaped (1 x N) expressions - a matrix row is a plain vector in the
    /// uBLAS semantics this type preserves - are accepted transposed.
    template<class TDerived>
    requires (Internals::IsRowShapedEigenExpression<TDerived, N>)
    array_1d(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.transpose().template cast<T>()) {}

    /// Construction from any uBLAS vector expression (interoperability with
    /// the dynamic uBLAS Vector and the uBLAS proxies/expressions).
    template<class TExpression>
    array_1d(const boost::numeric::ublas::vector_expression<TExpression>& rExpression) : BaseType()
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != N) << "Wrong size in the construction from a vector expression [ expression size = " << r_expression.size() << ", array size = " << N << " ]." << std::endl;
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] = r_expression(i);
        }
    }

    ///@}
    ///@name Operators
    ///@{

    // Element access (both operators, as in the uBLAS implementation)
    const_reference operator()(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    reference operator()(size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    const_reference operator[](size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }
    reference operator[](size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return BaseType::data()[i];
    }

    // Assignment
    array_1d& operator=(const array_1d& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// The assignment evaluates through a temporary (.eval()), reproducing
    /// the alias-safe uBLAS operator= semantics (e.g. x = trans-like
    /// expressions of x itself); the temporary-free fast path remains
    /// noalias(x) = expr, exactly as in uBLAS.
    template<class TDerived>
    requires (Internals::IsColumnShapedEigenExpression<TDerived, N> || Internals::IsRowShapedEigenExpression<TDerived, N>)
    array_1d& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        if constexpr (Internals::IsRowShapedEigenExpression<TDerived, N>) {
            BaseType::operator=(rExpression.transpose().template cast<T>().eval());
        } else {
            BaseType::operator=(rExpression.template cast<T>().eval());
        }
        return *this;
    }

    template<class TDerived>
    requires (Internals::IsColumnShapedEigenExpression<TDerived, N> || Internals::IsRowShapedEigenExpression<TDerived, N>)
    array_1d& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        if constexpr (Internals::IsRowShapedEigenExpression<TDerived, N>) {
            BaseType::operator+=(rExpression.transpose().template cast<T>().eval());
        } else {
            BaseType::operator+=(rExpression.template cast<T>().eval());
        }
        return *this;
    }

    template<class TDerived>
    requires (Internals::IsColumnShapedEigenExpression<TDerived, N> || Internals::IsRowShapedEigenExpression<TDerived, N>)
    array_1d& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        if constexpr (Internals::IsRowShapedEigenExpression<TDerived, N>) {
            BaseType::operator-=(rExpression.transpose().template cast<T>().eval());
        } else {
            BaseType::operator-=(rExpression.template cast<T>().eval());
        }
        return *this;
    }

    template<class TExpression>
    array_1d& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        // a uBLAS expression cannot reference this (Eigen) object, so no
        // protective temporary is needed
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != N) << "Wrong size in the assignment from a vector expression [ expression size = " << r_expression.size() << ", array size = " << N << " ]." << std::endl;
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] = r_expression(i);
        }
        return *this;
    }

    template<class TExpression>
    array_1d& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != N) << "Wrong size in the addition of a vector expression [ expression size = " << r_expression.size() << ", array size = " << N << " ]." << std::endl;
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] += r_expression(i);
        }
        return *this;
    }

    template<class TExpression>
    array_1d& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_expression = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_expression.size() != N) << "Wrong size in the subtraction of a vector expression [ expression size = " << r_expression.size() << ", array size = " << N << " ]." << std::endl;
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] -= r_expression(i);
        }
        return *this;
    }

    // Scalar operators keep the base class implementations (no aliasing there)
    using BaseType::operator*=;
    using BaseType::operator/=;

    /// Implicit conversion to the dynamic uBLAS vector, mirroring the
    /// implicit vector_expression conversion the uBLAS implementation offers
    /// (e.g. Vector v = some_array; or passing an array_1d to a Vector
    /// parameter).
    operator boost::numeric::ublas::vector<T>() const
    {
        boost::numeric::ublas::vector<T> result(N);
        for (size_type i = 0; i < N; ++i) {
            result[i] = (*this)[i];
        }
        return result;
    }

    /**
     * @brief Compares whether this array_1d is equal to the given array_1d.
     * @param rOther the array_1d to compare to
     * @return true if the two arrays are equal, false otherwise
     */
    bool operator==(const array_1d& rOther) const
    {
        return std::equal(begin(), end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    /// Interface-compatibility resize: the size is fixed, but !Preserve value-initializes.
    void resize(size_type Size, bool Preserve = true)
    {
        if (!Preserve) {
            std::fill(begin(), end(), value_type());
        }
    }

    array_1d& assign_temporary(array_1d& rOther)
    {
        swap(rOther);
        return *this;
    }

    // uBLAS assignment protocol (kept for generic code written against it)
    template<class TExpression>
    array_1d& assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this = rExpression);
    }
    template<class TExpression>
    array_1d& plus_assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this += rExpression);
    }
    template<class TExpression>
    array_1d& minus_assign(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        return (*this -= rExpression);
    }

    void swap(array_1d& rOther)
    {
        if (this != &rOther) {
            std::swap_ranges(begin(), end(), rOther.begin());
        }
    }
    friend void swap(array_1d& rFirst, array_1d& rSecond)
    {
        rFirst.swap(rSecond);
    }

    reference insert_element(size_type i, const_reference rValue)
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return ((*this)[i] = rValue);
    }
    void erase_element(size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i >= N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        (*this)[i] = value_type();
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
        return N;
    }

    /// Raw contiguous storage pointer (Eigen semantics: T* rather than the
    /// std::array& of the uBLAS implementation).
    using BaseType::data;

    iterator begin() { return BaseType::data(); }
    const_iterator begin() const { return BaseType::data(); }
    const_iterator cbegin() const { return BaseType::data(); }
    iterator end() { return BaseType::data() + N; }
    const_iterator end() const { return BaseType::data() + N; }
    const_iterator cend() const { return BaseType::data() + N; }

    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    ///@}
}; // Class array_1d

///@}
///@name Input and output
///@{

/// output stream function: reproduces the uBLAS text format [N](v0,v1,...)
template<class T, std::size_t N>
inline std::ostream& operator<<(std::ostream& rOStream, const array_1d<T, N>& rArray)
{
    rOStream << '[' << N << "](";
    for (std::size_t i = 0; i < N; ++i) {
        if (i > 0) {
            rOStream << ',';
        }
        rOStream << rArray[i];
    }
    rOStream << ')';
    return rOStream;
}

///@}

// The layout contract: Variable component access (variable.h), the
// DataTypeTraits contiguity used by the MPI buffers and the Point/Node
// coordinates all treat an array_1d<double, N> as a plain double[N] starting
// at the object address. Eigen::DontAlign in the base type is what disables
// the over-alignment Eigen would otherwise apply to the 4- and 6-sized
// arrays.
static_assert(sizeof(array_1d<double, 3>) == 3 * sizeof(double), "array_1d<double, 3> must be layout-compatible with double[3]");
static_assert(sizeof(array_1d<double, 4>) == 4 * sizeof(double), "array_1d<double, 4> must be layout-compatible with double[4]");
static_assert(sizeof(array_1d<double, 6>) == 6 * sizeof(double), "array_1d<double, 6> must be layout-compatible with double[6]");
static_assert(sizeof(array_1d<double, 9>) == 9 * sizeof(double), "array_1d<double, 9> must be layout-compatible with double[9]");
static_assert(alignof(array_1d<double, 4>) == alignof(double), "array_1d must not be over-aligned");

}  // namespace Kratos.

#else // uBLAS backend

// System includes
#include <string>
#include <iostream>
#include <array>
#include <algorithm>
#include <initializer_list>

// External	includes
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/detail/vector_assign.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp> // For container_const_reference etc.
#include <boost/numeric/ublas/functional.hpp> // For scalar_assign etc.

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type	Definitions
///@{

///@}
///@name	Enum's
///@{

///@}
///@name	Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short	class definition.
/** Detail class definition.
*/
template<class T,	std::size_t	N>
class	array_1d	: public boost::numeric::ublas::vector_expression< array_1d<T, N> >
{
public:
    ///@name Type	Definitions
    ///@{

    /// Pointer definition of	array_1d
    KRATOS_CLASS_POINTER_DEFINITION(array_1d);

    typedef std::size_t size_type;
    typedef	std::ptrdiff_t difference_type;
    typedef	T value_type;
    typedef	typename boost::numeric::ublas::type_traits<T>::const_reference const_reference;
    typedef	T &reference;
    typedef	std::array<T,N> array_type;
    typedef	T *pointer;
    typedef	array_1d<T, N> self_type;
    typedef	const boost::numeric::ublas::vector_reference<const self_type>	const_closure_type;
    typedef	boost::numeric::ublas::vector_reference<self_type>	closure_type;
    typedef	self_type vector_temporary_type;
    typedef	boost::numeric::ublas::dense_tag storage_category;

    ///@}
    ///@name Life	Cycle
    ///@{

    /// Default constructor.
    BOOST_UBLAS_INLINE
    array_1d ():
        boost::numeric::ublas::vector_expression<self_type> ()
    {
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size):
        boost::numeric::ublas::vector_expression<self_type> ()
    {
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size, value_type v):
        boost::numeric::ublas::vector_expression<self_type> ()
    {
        KRATOS_DEBUG_ERROR_IF(array_size>N) << "Given size is greater than the size of the array!" << std::endl;
        std::fill(data().begin(), data().begin() + array_size, v);
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (const std::initializer_list<value_type>& rInitList):
        boost::numeric::ublas::vector_expression<self_type> ()
    {
        KRATOS_DEBUG_ERROR_IF(rInitList.size()>N) << "Size of list greater than the size of the array!" << std::endl;
        std::copy(rInitList.begin(), rInitList.end(), data().begin());
    }

    BOOST_UBLAS_INLINE
    array_1d (size_type array_size,	const array_type & rdata):
        boost::numeric::ublas::vector_expression<self_type> (),
        data_ (rdata) {}

    BOOST_UBLAS_INLINE
    array_1d (const array_1d &v):
        boost::numeric::ublas::vector_expression<self_type> (),
        data_ (v.data_)	{}

    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d (const boost::numeric::ublas::vector_expression<AE> &ae) :
      boost::numeric::ublas::vector_expression<self_type> () // Initialize base class
    {
        boost::numeric::ublas::vector_assign<boost::numeric::ublas::scalar_assign> (*this, ae);
    }

    ///@}
    ///@name Operators
    ///@{

    // Element access
    BOOST_UBLAS_INLINE
    const_reference	operator ()	(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator () (size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }

    BOOST_UBLAS_INLINE
    const_reference	operator []	(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator [] (size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }

    // Assignment
    BOOST_UBLAS_INLINE
    array_1d &operator = (const array_1d &v)
    {
        data_ 	= v.data_;
        return *this;
    }

    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator = (const boost::numeric::ublas::vector_expression<AE>	&ae)
    {
        return assign (self_type	(ae));
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator +=	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        return assign (self_type	(*this + ae));
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator -=	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        return assign (self_type	(*this - ae));
    }
    template<class AT>
    BOOST_UBLAS_INLINE
    array_1d &operator /=	(const AT &at)
    {
        boost::numeric::ublas::vector_assign_scalar<boost::numeric::ublas::scalar_divides_assign> (*this, at);
        return *this;
    }

    /**
     * @brief Compares whether this array_1d is equal to the given array_1d.
     * @param v the array_1d to compare to
     * @return true if the two arrays are equal, false otherwise
     */
    BOOST_UBLAS_INLINE
    bool operator == (const array_1d &v) const
    {
        return std::equal (data_.begin(), data_.end(), v.data_.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    // Resizing
    BOOST_UBLAS_INLINE
    void resize	(size_type array_size, bool preserve = true)
    {
        if (!preserve)
            std::fill (data_.begin(), data_.end(), value_type	());
    }

    BOOST_UBLAS_INLINE
    array_1d &assign_temporary (array_1d &v)
    {
        swap (v);
        return *this;
    }

    template<class AT>
    BOOST_UBLAS_INLINE
    array_1d &operator *=	(const AT &at)
    {
        boost::numeric::ublas::vector_assign_scalar<boost::numeric::ublas::scalar_multiplies_assign> (*this, at);
        return *this;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &plus_assign	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        boost::numeric::ublas::vector_assign<boost::numeric::ublas::scalar_plus_assign> (*this, ae);
        return *this;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &assign (const boost::numeric::ublas::vector_expression<AE>	&ae)
    {
        boost::numeric::ublas::vector_assign<boost::numeric::ublas::scalar_assign> (*this, ae);
        return *this;
    }
    // Swapping
    BOOST_UBLAS_INLINE
    void swap (array_1d &v)
    {
        if (this !=	&v)
        {
            data ().swap (v.data ());
        }
    }
#ifndef	BOOST_UBLAS_NO_MEMBER_FRIENDS
    BOOST_UBLAS_INLINE
    friend void	swap (array_1d &v1, array_1d &v2)
    {
        v1.swap	(v2);
    }
#endif

    BOOST_UBLAS_INLINE
    reference insert_element (size_type i, const_reference t)
    {
        BOOST_UBLAS_CHECK (i < N, bad_index ());
        return (data_ [i] = t);
    }
    BOOST_UBLAS_INLINE
    void erase_element (size_type i)
    {
        BOOST_UBLAS_CHECK (i < N, bad_index ());
        data_ [i] = value_type/*zero*/();
    }
    BOOST_UBLAS_INLINE
    void clear ()
    {
        std::fill (data	().begin (), data ().end (), value_type	(0));
    }

    ///@}
    ///@name Access
    ///@{

private:
    typedef	typename array_type::const_iterator const_iterator_type;
    typedef	typename array_type::iterator iterator_type;

public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
    typedef	boost::numeric::ublas::indexed_iterator<self_type,	boost::numeric::ublas::dense_random_access_iterator_tag> iterator;
    typedef	boost::numeric::ublas::indexed_const_iterator<self_type, boost::numeric::ublas::dense_random_access_iterator_tag>	const_iterator;
#else
    class const_iterator;
    class iterator;
#endif

    // Element lookup
    BOOST_UBLAS_INLINE
    const_iterator find	(size_type i) const
    {
#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
        return const_iterator (*this, data ().begin	() + i);
#else
        return const_iterator (*this, i);
#endif
    }
    BOOST_UBLAS_INLINE
    iterator find (size_type i)
    {
#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
        return iterator	(*this,	data ().begin () + i);
#else
        return iterator	(*this,	i);
#endif
    }
    BOOST_UBLAS_INLINE
    size_type size () const
    {
        return N;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &minus_assign (const	boost::numeric::ublas::vector_expression<AE> &ae)
    {
        boost::numeric::ublas::vector_assign<boost::numeric::ublas::scalar_minus_assign>(*this,ae);
        return *this;
    }
    BOOST_UBLAS_INLINE
    const array_type &data () const
    {
        return data_;
    }
    BOOST_UBLAS_INLINE
    array_type &data ()
    {
        return data_;
    }

#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
    class const_iterator:
        public boost::numeric::ublas::container_const_reference<self_type>,
        public boost::numeric::ublas::random_access_iterator_base<typename boost::numeric::ublas::dense_random_access_iterator_tag, //dense_random_access_iterator_tag
        const_iterator, value_type, difference_type>
    {
    public:
        typedef	typename boost::numeric::ublas::dense_random_access_iterator_tag iterator_category; //dense_random_access_iterator_tag
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef	const_reference	reference;
#else
        typedef	typename array_1d::difference_type difference_type;
        typedef	typename array_1d::value_type	value_type;
        typedef	typename array_1d::const_reference reference;
        typedef	const typename array_1d::pointer pointer;
#endif

        BOOST_UBLAS_INLINE
        const_iterator ():
            boost::numeric::ublas::container_const_reference<self_type> (), it_ ()	{}
        BOOST_UBLAS_INLINE
        const_iterator (const self_type	&v,	const const_iterator_type &it):
            boost::numeric::ublas::container_const_reference<self_type> (v), it_ (it) {}
        BOOST_UBLAS_INLINE
#ifndef	BOOST_UBLAS_QUALIFIED_TYPENAME
        const_iterator (const iterator &it):
#else
        const_iterator (const typename self_type::iterator &it):
#endif
            boost::numeric::ublas::container_const_reference<self_type> (it ()), it_ (it.it_) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator &operator ++	()
        {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator --	()
        {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator +=	(difference_type n)
        {
            it_	+= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator -=	(difference_type n)
        {
            it_	-= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type	operator - (const const_iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference	operator * () const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index	());
            return *it_;
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index	() const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index	());
            return it_ - (*this) ().begin ().it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator &operator = (const const_iterator &it)
        {
            container_const_reference<self_type>::assign (&it ());
            it_	= it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const	const_iterator &it)	const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator <	(const const_iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        const_iterator_type	it_;

        friend class iterator;
    };
#endif

#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
    class iterator:
        public boost::numeric::ublas::container_reference<self_type>,
        public boost::numeric::ublas::random_access_iterator_base<typename boost::numeric::ublas::dense_random_access_iterator_tag, //dense_random_access_iterator_tag
        iterator, value_type, difference_type>
    {
    public:
        typedef	typename boost::numeric::ublas::dense_random_access_iterator_tag iterator_category; //dense_random_access_iterator_tag
#ifndef	BOOST_MSVC_STD_ITERATOR
        typedef	typename array_1d::difference_type difference_type;
        typedef	typename array_1d::value_type	value_type;
        typedef	typename array_1d::reference reference;
        typedef	typename array_1d::pointer pointer;
#endif

        BOOST_UBLAS_INLINE
        iterator ():
            boost::numeric::ublas::container_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        iterator (self_type	&v,	const iterator_type	&it):
            boost::numeric::ublas::container_reference<self_type> (v),	it_	(it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        iterator &operator ++ ()
        {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator -- ()
        {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator += (difference_type n)
        {
            it_	+= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator -= (difference_type n)
        {
            it_	-= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type	operator - (const iterator &it)	const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * ()	const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
            return *it_;
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index	() const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
            return it_ - (*this) ().begin ().it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        iterator &operator = (const	iterator &it)
        {
            container_reference<self_type>::assign (&it	());
            it_	= it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const	iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator <	(const iterator	&it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        iterator_type it_;
        friend class const_iterator;
    };
#endif

    BOOST_UBLAS_INLINE
    const_iterator begin ()	const
    {
        return find	(0);
    }
    BOOST_UBLAS_INLINE
    const_iterator end () const
    {
        return find	(data_.size	());
    }

    BOOST_UBLAS_INLINE
    iterator begin ()
    {
        return find	(0);
    }
    BOOST_UBLAS_INLINE
    iterator end ()
    {
        return find	(data_.size	());
    }

#ifdef BOOST_MSVC_STD_ITERATOR
    typedef	boost::numeric::ublas::reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
    typedef	boost::numeric::ublas::reverse_iterator_base<const_iterator> const_reverse_iterator;
#endif

    BOOST_UBLAS_INLINE
    const_reverse_iterator rbegin () const
    {
        return const_reverse_iterator (end ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator rend	() const
    {
        return const_reverse_iterator (begin ());
    }

#ifdef BOOST_MSVC_STD_ITERATOR
    typedef	boost::numeric::ublas::reverse_iterator_base<iterator,	value_type,	reference> reverse_iterator;
#else
    typedef	boost::numeric::ublas::reverse_iterator_base<iterator>	reverse_iterator;
#endif

    BOOST_UBLAS_INLINE
    reverse_iterator rbegin	()
    {
        return reverse_iterator	(end ());
    }
    BOOST_UBLAS_INLINE
    reverse_iterator rend ()
    {
        return reverse_iterator	(begin ());
    }
    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static	Member Variables
    ///@{

    ///@}
    ///@name Protected member	Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    array_type data_;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private	Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class	array_1d

///@}
///@name Type	Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace	Kratos.

#endif // KRATOS_USE_EIGEN_BACKEND

namespace AuxiliaryHashCombine
{
    /**
     * @brief This method creates an "unique" hash for the input value
     * @details It comes from boost, taken from here: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
     * @tparam TClassType The type of class to be hashed
     * @param rSeed This is the seed used to create the hash
     * @param rValue This is the value to be hashed
     * @todo Once the hashers and comparors are moved from key_hash.h, include key_hash and remove this. Right now there is a cross inclussion
     */
    template <class TClassType>
    inline void HashCombine(
        std::size_t& rSeed,
        const TClassType& rValue
        )
    {
        std::hash<TClassType> hasher;
        rSeed ^= hasher(rValue) + 0x9e3779b9 + (rSeed<<6) + (rSeed>>2);
    }
} /// namespace AuxiliaryHashCombine

namespace std
{
template<class T, std::size_t N>
struct hash<Kratos::array_1d<T,N>>
{
    std::size_t operator()(const Kratos::array_1d<T,N>& rArray) {
            std::size_t seed = 0;
            for (auto component : rArray) {AuxiliaryHashCombine::HashCombine(seed, component);}
            return seed;
        }
};
} /// namespace std
