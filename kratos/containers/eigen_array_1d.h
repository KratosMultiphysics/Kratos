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

// The Eigen implementation of array_1d (KRATOS_LINEAR_ALGEBRA_BACKEND=eigen,
// i.e. KRATOS_USE_EIGEN_BACKEND). This header is the Eigen arm of
// containers/array_1d.h and is not meant to be included directly.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <type_traits>

// External includes
#include <Eigen/Core>

// Project includes
#include "includes/define.h"
#include "includes/default_interface.h" // the Eigen-backed types and the uBLAS-style operations

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class array_1d
 * @brief Eigen-backed implementation of the Kratos fixed-size dense vector.
 * @details Selected by the Eigen linear-algebra backend
 * (KRATOS_USE_EIGEN_BACKEND); the uBLAS backend uses the historical
 * boost::numeric::ublas vector_expression implementation instead
 * (containers/ublas_array_1d.h).
 * It is a thin layer over a fixed-size Eigen column vector (every Eigen
 * expression is available on the type) that only adds what Eigen and the
 * Kratos Eigen plugins (includes/eigen_*_plugin.h) do not provide:
 * - Plain T[N] memory layout starting at the object address (no vtable, no
 *   padding, no over-alignment). Variable component access, DataTypeTraits
 *   contiguity and the MPI buffers rely on it, hence Eigen::DontAlign (pinned
 *   down by DenseBackendArray1dLayout in test_dense_backend_contract.cpp).
 * - The uBLAS-style constructors (size, size + value, initializer list),
 *   size_t-valued constexpr size(), checked operator()/operator[] (Debug),
 *   the data() storage view and a resize() that value-initializes when the
 *   content is not preserved.
 * - Element-wise scalar conversion on construction/assignment from other
 *   scalar types (Eigen requires an explicit cast).
 * - The uBLAS text format for streaming: [N](v0,v1,...).
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
        std::fill(this->begin(), this->begin() + Size, rValue);
    }

    explicit array_1d(const std::initializer_list<value_type>& rInitList) : BaseType()
    {
        KRATOS_DEBUG_ERROR_IF(rInitList.size() > N) << "Size of list greater than the size of the array!" << std::endl;
        std::copy(rInitList.begin(), rInitList.end(), this->begin());
    }

    array_1d(size_type Size, const array_type& rData) : BaseType()
    {
        std::copy(rData.begin(), rData.end(), this->begin());
    }

    array_1d(const array_1d& rOther) = default;

    /// Construction from any Eigen expression of compatible static size (the
    /// constraint keeps incompatible constructions out of overload sets). The
    /// cast mirrors the uBLAS element-wise conversion semantics for
    /// expressions of a different scalar type (e.g. array_1d<int, 3> from an
    /// array_1d<std::size_t, 3>).
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, N, 1>)
    array_1d(const Eigen::MatrixBase<TDerived>& rExpression) : BaseType(rExpression.template cast<T>()) {}

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

    array_1d& operator=(const array_1d& rOther) = default;

    /// Assignment from any compatible Eigen expression, converting the scalar
    /// type element-wise. Eigen evaluates through a temporary where aliasing
    /// matters (e.g. products); noalias(x) = expr skips it, exactly as in uBLAS.
    template<class TDerived>
    requires (Internals::IsBoundedCompatible<TDerived, N, 1>)
    array_1d& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        BaseType::operator=(rExpression.template cast<T>());
        return *this;
    }

    // The remaining operators (+=, -=, *=, /=) are Eigen's

    ///@}
    ///@name Operations
    ///@{

    /// Interface-compatibility resize: the size is fixed, but !Preserve value-initializes.
    void resize(size_type Size, bool Preserve = true)
    {
        if (!Preserve) {
            std::fill(this->begin(), this->end(), value_type());
        }
    }

    ///@}
    ///@name Access
    ///@{

    constexpr size_type size() const
    {
        return N;
    }

    /// Contiguous storage view (Kratos::span), with the uBLAS data() surface
    /// (begin/end as raw pointers, size, operator[]).
    Kratos::span<T> data() { return {BaseType::data(), N}; }
    Kratos::span<const T> data() const { return {BaseType::data(), N}; }

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

}  // namespace Kratos.
