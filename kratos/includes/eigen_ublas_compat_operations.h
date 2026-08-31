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

// External includes
#include <Eigen/Core>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Project includes
#include "includes/define.h"
#include "includes/eigen_compat_operations.h"

namespace Kratos {

// Forward declarations of the Eigen-backed fixed-size types (this header is
// pulled in by ublas_interface.h before containers/array_1d.h is complete).
template<class T, std::size_t N> class array_1d;
template<class T, std::size_t N> class EigenBoundedVector;
template<class T, std::size_t N1, std::size_t N2> class EigenBoundedMatrix;

///@name Mixed uBLAS/Eigen compatibility operations
///@{
// Under the Eigen backend the fixed-size dense types (array_1d,
// BoundedVector, BoundedMatrix) are Eigen-derived while the dynamic dense
// types (Vector, Matrix) remain uBLAS, so the ublas idioms mixing both worlds
// (noalias(rRHS) = prod(BoundedMatrix, array_1d), inner_prod(array_1d,
// Vector), Vector v = prod(M, a), ...) need dedicated overloads. The
// overload-resolution scheme follows eigen_compat_operations.h:
// - Functions taking expressions deduce one operand on Eigen::MatrixBase and
//   the other on ublas::vector_expression/matrix_expression, a combination
//   neither library provides, so there is never an ambiguity with the
//   pure-Eigen or pure-uBLAS overloads.
// - Mixed prod/outer_prod return an EAGER uBLAS dense result: every consumer
//   can absorb it (uBLAS targets natively, the Eigen-backed types through
//   their interop constructors/assignments), including the widespread
//   `Vector v = prod(...)` / `Matrix m = prod(...)` copy-initializations,
//   which could not be built from an Eigen expression. Mixed operator+/- keep
//   the Eigen side's plain type (fixed-size, stack allocated).
// - noalias() overloads on the CONCRETE target types (they win over the fully
//   generic ublas template by partial ordering) and return a proxy accepting
//   =, += and -= from both expression worlds. For uBLAS dense targets with a
//   uBLAS RHS the proxy reproduces the exact ublas::noalias semantics
//   (assign/plus_assign/minus_assign), so pre-existing all-uBLAS code is
//   unaffected; Eigen RHS are written through a zero-copy Eigen::Map over the
//   contiguous uBLAS storage.
///@}

namespace Internals {

/// noalias target adapter for the Eigen-derived dense types.
template<class TTarget>
class EigenTargetNoAliasProxy
{
public:
    explicit EigenTargetNoAliasProxy(TTarget& rTarget) : mrTarget(rTarget) {}

    template<class TDerived>
    TTarget& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() = rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TTarget& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() += rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TTarget& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        mrTarget.noalias() -= rExpression.derived();
        return mrTarget;
    }

    // uBLAS RHS: element-wise (a uBLAS expression cannot alias an Eigen target)
    template<class TExpression>
    TTarget& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size() != mrTarget.size()) << "noalias: vector expression of size " << r_e.size() << " assigned to a target of size " << mrTarget.size() << std::endl;
        for (std::size_t i = 0; i < mrTarget.size(); ++i) mrTarget[i] = r_e(i);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size() != mrTarget.size()) << "noalias: vector expression of size " << r_e.size() << " added to a target of size " << mrTarget.size() << std::endl;
        for (std::size_t i = 0; i < mrTarget.size(); ++i) mrTarget[i] += r_e(i);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size() != mrTarget.size()) << "noalias: vector expression of size " << r_e.size() << " subtracted from a target of size " << mrTarget.size() << std::endl;
        for (std::size_t i = 0; i < mrTarget.size(); ++i) mrTarget[i] -= r_e(i);
        return mrTarget;
    }

    template<class TExpression>
    TTarget& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size1() != mrTarget.size1() || r_e.size2() != mrTarget.size2()) << "noalias: matrix expression of size (" << r_e.size1() << "," << r_e.size2() << ") assigned to a target of size (" << mrTarget.size1() << "," << mrTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mrTarget.size1(); ++i)
            for (std::size_t j = 0; j < mrTarget.size2(); ++j)
                mrTarget(i, j) = r_e(i, j);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size1() != mrTarget.size1() || r_e.size2() != mrTarget.size2()) << "noalias: matrix expression of size (" << r_e.size1() << "," << r_e.size2() << ") added to a target of size (" << mrTarget.size1() << "," << mrTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mrTarget.size1(); ++i)
            for (std::size_t j = 0; j < mrTarget.size2(); ++j)
                mrTarget(i, j) += r_e(i, j);
        return mrTarget;
    }
    template<class TExpression>
    TTarget& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        const auto& r_e = rExpression();
        KRATOS_DEBUG_ERROR_IF(r_e.size1() != mrTarget.size1() || r_e.size2() != mrTarget.size2()) << "noalias: matrix expression of size (" << r_e.size1() << "," << r_e.size2() << ") subtracted from a target of size (" << mrTarget.size1() << "," << mrTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mrTarget.size1(); ++i)
            for (std::size_t j = 0; j < mrTarget.size2(); ++j)
                mrTarget(i, j) -= r_e(i, j);
        return mrTarget;
    }

private:
    TTarget& mrTarget;
};

/// noalias target adapter for the uBLAS dense vector: Eigen RHS through a
/// zero-copy Map, uBLAS RHS with the exact ublas::noalias semantics.
template<class TDataType>
class UblasVectorNoAliasProxy
{
public:
    using TargetType = boost::numeric::ublas::vector<TDataType>;
    using MapType = Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, 1>, Eigen::Unaligned>;

    explicit UblasVectorNoAliasProxy(TargetType& rTarget) : mrTarget(rTarget) {}

    template<class TDerived>
    TargetType& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size()).noalias() = rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TargetType& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size()).noalias() += rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TargetType& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size()).noalias() -= rExpression.derived();
        return mrTarget;
    }

    template<class TExpression>
    TargetType& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mrTarget.assign(rExpression);
        return mrTarget;
    }
    template<class TExpression>
    TargetType& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mrTarget.plus_assign(rExpression);
        return mrTarget;
    }
    template<class TExpression>
    TargetType& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mrTarget.minus_assign(rExpression);
        return mrTarget;
    }

private:
    TargetType& mrTarget;
};

/// noalias target adapter for the uBLAS dense (row-major) matrix.
template<class TDataType>
class UblasMatrixNoAliasProxy
{
public:
    using TargetType = boost::numeric::ublas::matrix<TDataType>;
    using MapType = Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned>;

    explicit UblasMatrixNoAliasProxy(TargetType& rTarget) : mrTarget(rTarget) {}

    template<class TDerived>
    TargetType& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size1(), mrTarget.size2()).noalias() = rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TargetType& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size1(), mrTarget.size2()).noalias() += rExpression.derived();
        return mrTarget;
    }
    template<class TDerived>
    TargetType& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        MapType(mrTarget.data().begin(), mrTarget.size1(), mrTarget.size2()).noalias() -= rExpression.derived();
        return mrTarget;
    }

    template<class TExpression>
    TargetType& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mrTarget.assign(rExpression);
        return mrTarget;
    }
    template<class TExpression>
    TargetType& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mrTarget.plus_assign(rExpression);
        return mrTarget;
    }
    template<class TExpression>
    TargetType& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mrTarget.minus_assign(rExpression);
        return mrTarget;
    }

private:
    TargetType& mrTarget;
};

/// noalias target adapter for the uBLAS vector-shaped PROXIES (matrix_row,
/// matrix_column, vector_range, ...), which are passed around by value: the
/// proxy view itself is copied, the referenced storage is written through.
/// uBLAS RHS reproduce the exact ublas::noalias semantics; Eigen RHS are
/// evaluated and copied element-wise.
template<class TProxy>
class UblasVectorProxyNoAliasProxy
{
public:
    explicit UblasVectorProxyNoAliasProxy(const TProxy& rTarget) : mTarget(rTarget) {}

    template<class TDerived>
    TProxy& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != mTarget.size()) << "noalias: expression of size " << evaluated.size() << " assigned to a proxy of size " << mTarget.size() << std::endl;
        for (std::size_t i = 0; i < mTarget.size(); ++i) mTarget(i) = evaluated(i);
        return mTarget;
    }
    template<class TDerived>
    TProxy& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != mTarget.size()) << "noalias: expression of size " << evaluated.size() << " added to a proxy of size " << mTarget.size() << std::endl;
        for (std::size_t i = 0; i < mTarget.size(); ++i) mTarget(i) += evaluated(i);
        return mTarget;
    }
    template<class TDerived>
    TProxy& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != mTarget.size()) << "noalias: expression of size " << evaluated.size() << " subtracted from a proxy of size " << mTarget.size() << std::endl;
        for (std::size_t i = 0; i < mTarget.size(); ++i) mTarget(i) -= evaluated(i);
        return mTarget;
    }

    template<class TExpression>
    TProxy& operator=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mTarget.assign(rExpression);
        return mTarget;
    }
    template<class TExpression>
    TProxy& operator+=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mTarget.plus_assign(rExpression);
        return mTarget;
    }
    template<class TExpression>
    TProxy& operator-=(const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
    {
        mTarget.minus_assign(rExpression);
        return mTarget;
    }

private:
    mutable TProxy mTarget;
};

/// Matrix-shaped counterpart (matrix_range, matrix_slice).
template<class TProxy>
class UblasMatrixProxyNoAliasProxy
{
public:
    explicit UblasMatrixProxyNoAliasProxy(const TProxy& rTarget) : mTarget(rTarget) {}

    template<class TDerived>
    TProxy& operator=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.rows()) != mTarget.size1() || static_cast<std::size_t>(evaluated.cols()) != mTarget.size2()) << "noalias: expression of size (" << evaluated.rows() << "," << evaluated.cols() << ") assigned to a proxy of size (" << mTarget.size1() << "," << mTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mTarget.size1(); ++i)
            for (std::size_t j = 0; j < mTarget.size2(); ++j)
                mTarget(i, j) = evaluated(i, j);
        return mTarget;
    }
    template<class TDerived>
    TProxy& operator+=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.rows()) != mTarget.size1() || static_cast<std::size_t>(evaluated.cols()) != mTarget.size2()) << "noalias: expression of size (" << evaluated.rows() << "," << evaluated.cols() << ") added to a proxy of size (" << mTarget.size1() << "," << mTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mTarget.size1(); ++i)
            for (std::size_t j = 0; j < mTarget.size2(); ++j)
                mTarget(i, j) += evaluated(i, j);
        return mTarget;
    }
    template<class TDerived>
    TProxy& operator-=(const Eigen::MatrixBase<TDerived>& rExpression)
    {
        const auto evaluated = rExpression.derived().eval();
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.rows()) != mTarget.size1() || static_cast<std::size_t>(evaluated.cols()) != mTarget.size2()) << "noalias: expression of size (" << evaluated.rows() << "," << evaluated.cols() << ") subtracted from a proxy of size (" << mTarget.size1() << "," << mTarget.size2() << ")" << std::endl;
        for (std::size_t i = 0; i < mTarget.size1(); ++i)
            for (std::size_t j = 0; j < mTarget.size2(); ++j)
                mTarget(i, j) -= evaluated(i, j);
        return mTarget;
    }

    template<class TExpression>
    TProxy& operator=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mTarget.assign(rExpression);
        return mTarget;
    }
    template<class TExpression>
    TProxy& operator+=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mTarget.plus_assign(rExpression);
        return mTarget;
    }
    template<class TExpression>
    TProxy& operator-=(const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
    {
        mTarget.minus_assign(rExpression);
        return mTarget;
    }

private:
    mutable TProxy mTarget;
};

/// Concrete contiguous uBLAS containers (the dynamic vector/matrix over an
/// unbounded_array): their storage can be viewed zero-copy through Eigen maps,
/// so the mixed products below can run on Eigen's vectorized kernels. uBLAS
/// proxy/lazy expressions do not qualify and take the scalar fallback paths.
template<class T> struct IsDenseUblasVector : std::false_type {};
template<class T> struct IsDenseUblasVector<boost::numeric::ublas::vector<T>> : std::true_type {};
template<class T> struct IsDenseUblasMatrix : std::false_type {};
template<class T> struct IsDenseUblasMatrix<boost::numeric::ublas::matrix<T>> : std::true_type {};

template<class T>
inline auto MapOfUblas(const boost::numeric::ublas::vector<T>& rV)
{
    using MapType = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Unaligned>;
    return MapType(rV.data().begin(), static_cast<Eigen::Index>(rV.size()));
}

template<class T>
inline auto MapOfUblas(const boost::numeric::ublas::matrix<T>& rM)
{
    using MapType = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned>;
    return MapType(rM.data().begin(), static_cast<Eigen::Index>(rM.size1()), static_cast<Eigen::Index>(rM.size2()));
}

template<class T>
inline auto MutableMapOfUblas(boost::numeric::ublas::vector<T>& rV)
{
    using MapType = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Unaligned>;
    return MapType(rV.data().begin(), static_cast<Eigen::Index>(rV.size()));
}

template<class T>
inline auto MutableMapOfUblas(boost::numeric::ublas::matrix<T>& rM)
{
    using MapType = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned>;
    return MapType(rM.data().begin(), static_cast<Eigen::Index>(rM.size1()), static_cast<Eigen::Index>(rM.size2()));
}

} // namespace Internals

///@name noalias
///@{

// The uBLAS proxies are returned by value, so the ublas noalias(const C&)
// overload is the one that binds; these more specialized overloads take over
// for the proxy types so an Eigen RHS also works.
template<class TMatrix>
inline auto noalias(const boost::numeric::ublas::matrix_row<TMatrix>& rTarget)
{
    return Internals::UblasVectorProxyNoAliasProxy<boost::numeric::ublas::matrix_row<TMatrix>>(rTarget);
}

template<class TMatrix>
inline auto noalias(const boost::numeric::ublas::matrix_column<TMatrix>& rTarget)
{
    return Internals::UblasVectorProxyNoAliasProxy<boost::numeric::ublas::matrix_column<TMatrix>>(rTarget);
}

template<class TVector>
inline auto noalias(const boost::numeric::ublas::vector_range<TVector>& rTarget)
{
    return Internals::UblasVectorProxyNoAliasProxy<boost::numeric::ublas::vector_range<TVector>>(rTarget);
}

template<class TMatrix>
inline auto noalias(const boost::numeric::ublas::matrix_range<TMatrix>& rTarget)
{
    return Internals::UblasMatrixProxyNoAliasProxy<boost::numeric::ublas::matrix_range<TMatrix>>(rTarget);
}

template<class T, std::size_t N>
inline auto noalias(array_1d<T, N>& rTarget)
{
    return Internals::EigenTargetNoAliasProxy<array_1d<T, N>>(rTarget);
}

template<class T, std::size_t N>
inline auto noalias(EigenBoundedVector<T, N>& rTarget)
{
    return Internals::EigenTargetNoAliasProxy<EigenBoundedVector<T, N>>(rTarget);
}

template<class T, std::size_t N1, std::size_t N2>
inline auto noalias(EigenBoundedMatrix<T, N1, N2>& rTarget)
{
    return Internals::EigenTargetNoAliasProxy<EigenBoundedMatrix<T, N1, N2>>(rTarget);
}

/// noalias on a writable Eigen expression proxy (e.g. row(rBoundedMatrix, i)
/// or subrange(rArray, ...) results, which are passed around by value).
/// Following the Eigen convention for write-through proxies the target is
/// taken by const reference and cast back; the requires-clause limits this to
/// lvalue (writable) Eigen expressions, and the concrete overloads above
/// remain preferred for the plain types.
template<class TDerived>
requires (static_cast<bool>(Eigen::internal::is_lvalue<TDerived>::value))
inline auto noalias(const Eigen::MatrixBase<TDerived>& rTarget)
{
    return Internals::EigenTargetNoAliasProxy<TDerived>(const_cast<Eigen::MatrixBase<TDerived>&>(rTarget).derived());
}

template<class T>
inline auto noalias(boost::numeric::ublas::vector<T>& rTarget)
{
    return Internals::UblasVectorNoAliasProxy<T>(rTarget);
}

template<class T>
inline auto noalias(boost::numeric::ublas::matrix<T>& rTarget)
{
    return Internals::UblasMatrixNoAliasProxy<T>(rTarget);
}

///@}
///@name Mixed products
///@{

/// (Eigen matrix) x (uBLAS vector expression) -> uBLAS dense vector.
/// The result stays an eager uBLAS container (see the design note above), but
/// when the uBLAS operand is a concrete dense vector the computation runs on
/// Eigen's vectorized product kernels over zero-copy maps.
template<class TDerived, class TExpression>
inline auto prod(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::vector_expression<TExpression>& rX)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& a = rA.derived();
    const auto& r_x = rX();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(a.cols()) != r_x.size()) << "prod: incompatible sizes [ matrix columns = " << a.cols() << ", vector size = " << r_x.size() << " ]." << std::endl;
    boost::numeric::ublas::vector<scalar_t> result(a.rows());
    if constexpr (Internals::IsDenseUblasVector<TExpression>::value) {
        Internals::MutableMapOfUblas(result).noalias() = a * Internals::MapOfUblas(r_x);
    } else {
        const auto a_eval = a.eval(); // coefficient loops must not re-evaluate an expression per access
        for (Eigen::Index i = 0; i < a_eval.rows(); ++i) {
            scalar_t aux = scalar_t();
            for (Eigen::Index k = 0; k < a_eval.cols(); ++k) {
                aux += a_eval(i, k) * r_x(k);
            }
            result[i] = aux;
        }
    }
    return result;
}

/// (uBLAS matrix expression) x (Eigen vector or matrix). A single overload
/// dispatching on the compile-time shape of the Eigen side (both shapes
/// deduce as Eigen::MatrixBase, so two separate templates would collide):
/// vector RHS -> uBLAS dense vector, matrix RHS -> uBLAS dense matrix.
template<class TExpression, class TDerived>
inline auto prod(const boost::numeric::ublas::matrix_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& r_a = rA();
    if constexpr (TDerived::ColsAtCompileTime == 1) {
        KRATOS_DEBUG_ERROR_IF(r_a.size2() != static_cast<std::size_t>(rB.size())) << "prod: incompatible sizes [ matrix columns = " << r_a.size2() << ", vector size = " << rB.size() << " ]." << std::endl;
        boost::numeric::ublas::vector<scalar_t> result(r_a.size1());
        if constexpr (Internals::IsDenseUblasMatrix<TExpression>::value) {
            Internals::MutableMapOfUblas(result).noalias() = Internals::MapOfUblas(r_a) * rB.derived();
        } else {
            const auto x = rB.derived().eval();
            for (std::size_t i = 0; i < r_a.size1(); ++i) {
                scalar_t aux = scalar_t();
                for (std::size_t k = 0; k < r_a.size2(); ++k) {
                    aux += r_a(i, k) * x[k];
                }
                result[i] = aux;
            }
        }
        return result;
    } else {
        KRATOS_DEBUG_ERROR_IF(r_a.size2() != static_cast<std::size_t>(rB.rows())) << "prod: incompatible sizes [ left columns = " << r_a.size2() << ", right rows = " << rB.rows() << " ]." << std::endl;
        boost::numeric::ublas::matrix<scalar_t> result(r_a.size1(), rB.cols());
        if constexpr (Internals::IsDenseUblasMatrix<TExpression>::value) {
            Internals::MutableMapOfUblas(result).noalias() = Internals::MapOfUblas(r_a) * rB.derived();
        } else {
            const auto b = rB.derived().eval();
            for (std::size_t i = 0; i < r_a.size1(); ++i) {
                for (Eigen::Index j = 0; j < b.cols(); ++j) {
                    scalar_t aux = scalar_t();
                    for (std::size_t k = 0; k < r_a.size2(); ++k) {
                        aux += r_a(i, k) * b(k, j);
                    }
                    result(i, j) = aux;
                }
            }
        }
        return result;
    }
}

/// (Eigen matrix or vector) x (uBLAS matrix expression). Dispatches on the
/// compile-time shape of the Eigen side: a (column) vector v gives the ublas
/// prod(v, M) semantics (v^T M -> uBLAS dense vector), a matrix gives the
/// matrix product (-> uBLAS dense matrix).
template<class TDerived, class TExpression>
inline auto prod(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::matrix_expression<TExpression>& rB)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& r_b = rB();
    if constexpr (TDerived::ColsAtCompileTime == 1) {
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(rA.size()) != r_b.size1()) << "prod: incompatible sizes [ vector size = " << rA.size() << ", matrix rows = " << r_b.size1() << " ]." << std::endl;
        boost::numeric::ublas::vector<scalar_t> result(r_b.size2());
        if constexpr (Internals::IsDenseUblasMatrix<TExpression>::value) {
            // ublas prod(v, M) semantics: result = M^T v
            Internals::MutableMapOfUblas(result).noalias() = Internals::MapOfUblas(r_b).transpose() * rA.derived();
        } else {
            const auto a = rA.derived().eval();
            for (std::size_t j = 0; j < r_b.size2(); ++j) {
                scalar_t aux = scalar_t();
                for (std::size_t k = 0; k < r_b.size1(); ++k) {
                    aux += a[k] * r_b(k, j);
                }
                result[j] = aux;
            }
        }
        return result;
    } else {
        KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(rA.cols()) != r_b.size1()) << "prod: incompatible sizes [ left columns = " << rA.cols() << ", right rows = " << r_b.size1() << " ]." << std::endl;
        boost::numeric::ublas::matrix<scalar_t> result(rA.rows(), r_b.size2());
        if constexpr (Internals::IsDenseUblasMatrix<TExpression>::value) {
            Internals::MutableMapOfUblas(result).noalias() = rA.derived() * Internals::MapOfUblas(r_b);
        } else {
            const auto a = rA.derived().eval();
            for (Eigen::Index i = 0; i < a.rows(); ++i) {
                for (std::size_t j = 0; j < r_b.size2(); ++j) {
                    scalar_t aux = scalar_t();
                    for (Eigen::Index k = 0; k < a.cols(); ++k) {
                        aux += a(i, k) * r_b(k, j);
                    }
                    result(i, j) = aux;
                }
            }
        }
        return result;
    }
}

/// Mixed scalar products
template<class TDerived, class TExpression>
inline auto inner_prod(const Eigen::MatrixBase<TDerived>& rX, const boost::numeric::ublas::vector_expression<TExpression>& rY)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& r_y = rY();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(rX.size()) != r_y.size()) << "inner_prod: incompatible sizes [ " << rX.size() << " and " << r_y.size() << " ]." << std::endl;
    if constexpr (Internals::IsDenseUblasVector<TExpression>::value) {
        return scalar_t(rX.derived().dot(Internals::MapOfUblas(r_y)));
    } else {
        const auto x = rX.derived().eval();
        scalar_t result = scalar_t();
        for (std::size_t i = 0; i < r_y.size(); ++i) {
            result += x[i] * r_y(i);
        }
        return result;
    }
}

template<class TExpression, class TDerived>
inline auto inner_prod(const boost::numeric::ublas::vector_expression<TExpression>& rX, const Eigen::MatrixBase<TDerived>& rY)
{
    return inner_prod(rY, rX);
}

/// Mixed outer products -> uBLAS dense matrix
template<class TDerived, class TExpression>
inline auto outer_prod(const Eigen::MatrixBase<TDerived>& rX, const boost::numeric::ublas::vector_expression<TExpression>& rY)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& r_y = rY();
    boost::numeric::ublas::matrix<scalar_t> result(rX.size(), r_y.size());
    // The map fast path needs a column-shaped Eigen operand; the scalar
    // fallback flat-indexes and also covers row-shaped views.
    if constexpr (Internals::IsDenseUblasVector<TExpression>::value && TDerived::ColsAtCompileTime == 1) {
        Internals::MutableMapOfUblas(result).noalias() = rX.derived() * Internals::MapOfUblas(r_y).transpose();
    } else {
        const auto x = rX.derived().eval();
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            for (std::size_t j = 0; j < r_y.size(); ++j) {
                result(i, j) = x[i] * r_y(j);
            }
        }
    }
    return result;
}

template<class TExpression, class TDerived>
inline auto outer_prod(const boost::numeric::ublas::vector_expression<TExpression>& rX, const Eigen::MatrixBase<TDerived>& rY)
{
    using scalar_t = typename TDerived::Scalar;
    const auto& r_x = rX();
    boost::numeric::ublas::matrix<scalar_t> result(r_x.size(), rY.size());
    // The map fast path needs a column-shaped Eigen operand; the scalar
    // fallback flat-indexes and also covers row-shaped views.
    if constexpr (Internals::IsDenseUblasVector<TExpression>::value && TDerived::ColsAtCompileTime == 1) {
        Internals::MutableMapOfUblas(result).noalias() = Internals::MapOfUblas(r_x) * rY.derived().transpose();
    } else {
        const auto y = rY.derived().eval();
        for (std::size_t i = 0; i < r_x.size(); ++i) {
            for (Eigen::Index j = 0; j < y.size(); ++j) {
                result(i, j) = r_x(i) * y[j];
            }
        }
    }
    return result;
}

///@}
///@name Mixed additive operations
///@{

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject operator+(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::vector_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.size()) != r_b.size()) << "operator+: incompatible sizes [ " << result.size() << " and " << r_b.size() << " ]." << std::endl;
    for (std::size_t i = 0; i < r_b.size(); ++i) {
        result[i] += r_b(i);
    }
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject operator+(const boost::numeric::ublas::vector_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    return rB + rA;
}

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject operator-(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::vector_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.size()) != r_b.size()) << "operator-: incompatible sizes [ " << result.size() << " and " << r_b.size() << " ]." << std::endl;
    for (std::size_t i = 0; i < r_b.size(); ++i) {
        result[i] -= r_b(i);
    }
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject operator-(const boost::numeric::ublas::vector_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    typename TDerived::PlainObject result = rB.derived();
    const auto& r_a = rA();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.size()) != r_a.size()) << "operator-: incompatible sizes [ " << result.size() << " and " << r_a.size() << " ]." << std::endl;
    for (std::size_t i = 0; i < r_a.size(); ++i) {
        result[i] = r_a(i) - result[i];
    }
    return result;
}

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject operator+(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::matrix_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.rows()) != r_b.size1() || static_cast<std::size_t>(result.cols()) != r_b.size2()) << "operator+: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_b.size1(); ++i)
        for (std::size_t j = 0; j < r_b.size2(); ++j)
            result(i, j) += r_b(i, j);
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject operator+(const boost::numeric::ublas::matrix_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    return rB + rA;
}

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject operator-(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::matrix_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.rows()) != r_b.size1() || static_cast<std::size_t>(result.cols()) != r_b.size2()) << "operator-: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_b.size1(); ++i)
        for (std::size_t j = 0; j < r_b.size2(); ++j)
            result(i, j) -= r_b(i, j);
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject operator-(const boost::numeric::ublas::matrix_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    typename TDerived::PlainObject result = rB.derived();
    const auto& r_a = rA();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(result.rows()) != r_a.size1() || static_cast<std::size_t>(result.cols()) != r_a.size2()) << "operator-: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_a.size1(); ++i)
        for (std::size_t j = 0; j < r_a.size2(); ++j)
            result(i, j) = r_a(i, j) - result(i, j);
    return result;
}

/// Compound assignment of Eigen expressions into the uBLAS dense containers
/// (e.g. rLeftHandSideMatrix += <bounded expression>), through a zero-copy Map.
template<class TDataType, class TDerived>
inline boost::numeric::ublas::vector<TDataType>& operator+=(boost::numeric::ublas::vector<TDataType>& rTarget, const Eigen::MatrixBase<TDerived>& rExpression)
{
    Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, 1>, Eigen::Unaligned>(rTarget.data().begin(), rTarget.size()) += rExpression.derived();
    return rTarget;
}

template<class TDataType, class TDerived>
inline boost::numeric::ublas::vector<TDataType>& operator-=(boost::numeric::ublas::vector<TDataType>& rTarget, const Eigen::MatrixBase<TDerived>& rExpression)
{
    Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, 1>, Eigen::Unaligned>(rTarget.data().begin(), rTarget.size()) -= rExpression.derived();
    return rTarget;
}

template<class TDataType, class TDerived>
inline boost::numeric::ublas::matrix<TDataType>& operator+=(boost::numeric::ublas::matrix<TDataType>& rTarget, const Eigen::MatrixBase<TDerived>& rExpression)
{
    Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned>(rTarget.data().begin(), rTarget.size1(), rTarget.size2()) += rExpression.derived();
    return rTarget;
}

template<class TDataType, class TDerived>
inline boost::numeric::ublas::matrix<TDataType>& operator-=(boost::numeric::ublas::matrix<TDataType>& rTarget, const Eigen::MatrixBase<TDerived>& rExpression)
{
    Eigen::Map<Eigen::Matrix<TDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned>(rTarget.data().begin(), rTarget.size1(), rTarget.size2()) -= rExpression.derived();
    return rTarget;
}

///@}
///@name Mixed element-wise operations
///@{

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject element_prod(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::vector_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    for (std::size_t i = 0; i < r_b.size(); ++i) result[i] *= r_b(i);
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject element_prod(const boost::numeric::ublas::vector_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    return element_prod(rB, rA);
}

template<class TDerived, class TExpression>
inline typename TDerived::PlainObject element_div(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::vector_expression<TExpression>& rB)
{
    typename TDerived::PlainObject result = rA.derived();
    const auto& r_b = rB();
    for (std::size_t i = 0; i < r_b.size(); ++i) result[i] /= r_b(i);
    return result;
}

template<class TExpression, class TDerived>
inline typename TDerived::PlainObject element_div(const boost::numeric::ublas::vector_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    typename TDerived::PlainObject result = rB.derived();
    const auto& r_a = rA();
    for (std::size_t i = 0; i < r_a.size(); ++i) result[i] = r_a(i) / result[i];
    return result;
}

/// Mixed matrix element-wise products -> uBLAS dense matrix (either side can
/// absorb it)
template<class TDerived, class TExpression>
inline auto element_prod(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::matrix_expression<TExpression>& rB)
{
    using scalar_t = typename TDerived::Scalar;
    const auto a = rA.derived().eval();
    const auto& r_b = rB();
    boost::numeric::ublas::matrix<scalar_t> result(r_b.size1(), r_b.size2());
    for (std::size_t i = 0; i < r_b.size1(); ++i)
        for (std::size_t j = 0; j < r_b.size2(); ++j)
            result(i, j) = a(i, j) * r_b(i, j);
    return result;
}

template<class TExpression, class TDerived>
inline auto element_prod(const boost::numeric::ublas::matrix_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    return element_prod(rB, rA);
}

template<class TDerived, class TExpression>
inline auto element_div(const Eigen::MatrixBase<TDerived>& rA, const boost::numeric::ublas::matrix_expression<TExpression>& rB)
{
    using scalar_t = typename TDerived::Scalar;
    const auto a = rA.derived().eval();
    const auto& r_b = rB();
    boost::numeric::ublas::matrix<scalar_t> result(r_b.size1(), r_b.size2());
    for (std::size_t i = 0; i < r_b.size1(); ++i)
        for (std::size_t j = 0; j < r_b.size2(); ++j)
            result(i, j) = a(i, j) / r_b(i, j);
    return result;
}

template<class TExpression, class TDerived>
inline auto element_div(const boost::numeric::ublas::matrix_expression<TExpression>& rA, const Eigen::MatrixBase<TDerived>& rB)
{
    using scalar_t = typename TDerived::Scalar;
    const auto b = rB.derived().eval();
    const auto& r_a = rA();
    boost::numeric::ublas::matrix<scalar_t> result(r_a.size1(), r_a.size2());
    for (std::size_t i = 0; i < r_a.size1(); ++i)
        for (std::size_t j = 0; j < r_a.size2(); ++j)
            result(i, j) = r_a(i, j) / b(i, j);
    return result;
}

///@}
///@name Compound assignment into uBLAS proxies and writable Eigen expressions
///@{

// The uBLAS proxies (project/subrange/row/column results) are passed around
// by value; taking them by value here lets the operators bind to those
// rvalues while still writing through the referenced storage.

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_range<TMatrix> operator+=(boost::numeric::ublas::matrix_range<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.rows()) != Target.size1() || static_cast<std::size_t>(evaluated.cols()) != Target.size2()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size1(); ++i)
        for (std::size_t j = 0; j < Target.size2(); ++j)
            Target(i, j) += evaluated(i, j);
    return Target;
}

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_range<TMatrix> operator-=(boost::numeric::ublas::matrix_range<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.rows()) != Target.size1() || static_cast<std::size_t>(evaluated.cols()) != Target.size2()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size1(); ++i)
        for (std::size_t j = 0; j < Target.size2(); ++j)
            Target(i, j) -= evaluated(i, j);
    return Target;
}

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_row<TMatrix> operator+=(boost::numeric::ublas::matrix_row<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) += evaluated(i);
    return Target;
}

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_row<TMatrix> operator-=(boost::numeric::ublas::matrix_row<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) -= evaluated(i);
    return Target;
}

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_column<TMatrix> operator+=(boost::numeric::ublas::matrix_column<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) += evaluated(i);
    return Target;
}

template<class TMatrix, class TDerived>
inline boost::numeric::ublas::matrix_column<TMatrix> operator-=(boost::numeric::ublas::matrix_column<TMatrix> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) -= evaluated(i);
    return Target;
}

template<class TVector, class TDerived>
inline boost::numeric::ublas::vector_range<TVector> operator+=(boost::numeric::ublas::vector_range<TVector> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) += evaluated(i);
    return Target;
}

template<class TVector, class TDerived>
inline boost::numeric::ublas::vector_range<TVector> operator-=(boost::numeric::ublas::vector_range<TVector> Target, const Eigen::MatrixBase<TDerived>& rExpression)
{
    const auto evaluated = rExpression.derived().eval();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(evaluated.size()) != Target.size()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < Target.size(); ++i)
        Target(i) -= evaluated(i);
    return Target;
}

// Compound assignment of a uBLAS expression into a writable Eigen expression
// (e.g. row(rBoundedMatrix, i) += <ublas row expression>). Following the
// Eigen convention for write-through proxies the target is taken by const
// reference and cast back; the requires-clause limits this to lvalue
// (writable) Eigen expressions.
template<class TDerived, class TExpression>
requires (static_cast<bool>(Eigen::internal::is_lvalue<TDerived>::value))
inline TDerived& operator+=(const Eigen::MatrixBase<TDerived>& rTarget, const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
{
    auto& r_target = const_cast<Eigen::MatrixBase<TDerived>&>(rTarget).derived();
    const auto& r_e = rExpression();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(r_target.size()) != r_e.size()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_e.size(); ++i)
        r_target(i) += r_e(i);
    return r_target;
}

template<class TDerived, class TExpression>
requires (static_cast<bool>(Eigen::internal::is_lvalue<TDerived>::value))
inline TDerived& operator-=(const Eigen::MatrixBase<TDerived>& rTarget, const boost::numeric::ublas::vector_expression<TExpression>& rExpression)
{
    auto& r_target = const_cast<Eigen::MatrixBase<TDerived>&>(rTarget).derived();
    const auto& r_e = rExpression();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(r_target.size()) != r_e.size()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_e.size(); ++i)
        r_target(i) -= r_e(i);
    return r_target;
}

template<class TDerived, class TExpression>
requires (static_cast<bool>(Eigen::internal::is_lvalue<TDerived>::value))
inline TDerived& operator+=(const Eigen::MatrixBase<TDerived>& rTarget, const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
{
    auto& r_target = const_cast<Eigen::MatrixBase<TDerived>&>(rTarget).derived();
    const auto& r_e = rExpression();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(r_target.rows()) != r_e.size1() || static_cast<std::size_t>(r_target.cols()) != r_e.size2()) << "operator+=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_e.size1(); ++i)
        for (std::size_t j = 0; j < r_e.size2(); ++j)
            r_target(i, j) += r_e(i, j);
    return r_target;
}

template<class TDerived, class TExpression>
requires (static_cast<bool>(Eigen::internal::is_lvalue<TDerived>::value))
inline TDerived& operator-=(const Eigen::MatrixBase<TDerived>& rTarget, const boost::numeric::ublas::matrix_expression<TExpression>& rExpression)
{
    auto& r_target = const_cast<Eigen::MatrixBase<TDerived>&>(rTarget).derived();
    const auto& r_e = rExpression();
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(r_target.rows()) != r_e.size1() || static_cast<std::size_t>(r_target.cols()) != r_e.size2()) << "operator-=: incompatible sizes." << std::endl;
    for (std::size_t i = 0; i < r_e.size1(); ++i)
        for (std::size_t j = 0; j < r_e.size2(); ++j)
            r_target(i, j) -= r_e(i, j);
    return r_target;
}

///@}
///@name Proxy accessors on the Eigen-backed bounded types
///@{

/// Row proxy: vector-shaped (transposed) so it reads into and assigns from
/// column-vector-shaped objects with the ublas row() semantics; writable.
template<class T, std::size_t N1, std::size_t N2>
inline auto row(EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
}

template<class T, std::size_t N1, std::size_t N2>
inline auto row(const EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t I)
{
    return rM.row(I).transpose();
}

/// Column proxy (readable and writable, as ublas column()).
template<class T, std::size_t N1, std::size_t N2>
inline auto column(EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t J)
{
    return rM.col(J);
}

template<class T, std::size_t N1, std::size_t N2>
inline auto column(const EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t J)
{
    return rM.col(J);
}

/// Vector subrange [Low, High) proxies (readable and writable).
template<class T, std::size_t N>
inline auto subrange(array_1d<T, N>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

template<class T, std::size_t N>
inline auto subrange(const array_1d<T, N>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

template<class T, std::size_t N>
inline auto subrange(EigenBoundedVector<T, N>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

template<class T, std::size_t N>
inline auto subrange(const EigenBoundedVector<T, N>& rV, const std::size_t Low, const std::size_t High)
{
    return rV.segment(Low, High - Low);
}

/// Matrix subrange [Row1, Row2) x [Col1, Col2) proxies (readable and writable).
template<class T, std::size_t N1, std::size_t N2>
inline auto subrange(EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

template<class T, std::size_t N1, std::size_t N2>
inline auto subrange(const EigenBoundedMatrix<T, N1, N2>& rM, const std::size_t Row1, const std::size_t Row2, const std::size_t Col1, const std::size_t Col2)
{
    return rM.block(Row1, Col1, Row2 - Row1, Col2 - Col1);
}

/// uBLAS trans() on a *vector* is the identity ((trans v)[i] = v[i]); see the
/// EigenVector overload in eigen_compat_operations.h. These cover the
/// fixed-size vector types (N x 1 *matrices* keep the real transpose).
template<class T, std::size_t N>
inline const EigenBoundedVector<T, N>& trans(const EigenBoundedVector<T, N>& rV)
{
    return rV;
}

template<class T, std::size_t N>
inline const array_1d<T, N>& trans(const array_1d<T, N>& rV)
{
    return rV;
}

/// ublas project() on the Eigen-backed types: a ublas range maps to a
/// segment/block view (readable and writable, exactly as the ublas proxies).
/// Ranges are preprocess()-ed against the target extent so the range::all()
/// sentinel resolves to the full extent, as in the ublas proxies.
template<class TDataType>
inline auto project(EigenVector<TDataType>& rV, const boost::numeric::ublas::range& rRange)
{
    const auto processed = rRange.preprocess(rV.size());
    return rV.segment(processed.start(), processed.size());
}

template<class TDataType>
inline auto project(const EigenVector<TDataType>& rV, const boost::numeric::ublas::range& rRange)
{
    const auto processed = rRange.preprocess(rV.size());
    return rV.segment(processed.start(), processed.size());
}

template<class TDataType>
inline auto project(EigenMatrix<TDataType>& rM, const boost::numeric::ublas::range& rRowRange, const boost::numeric::ublas::range& rColRange)
{
    const auto processed_rows = rRowRange.preprocess(rM.rows());
    const auto processed_cols = rColRange.preprocess(rM.cols());
    return rM.block(processed_rows.start(), processed_cols.start(), processed_rows.size(), processed_cols.size());
}

template<class TDataType>
inline auto project(const EigenMatrix<TDataType>& rM, const boost::numeric::ublas::range& rRowRange, const boost::numeric::ublas::range& rColRange)
{
    const auto processed_rows = rRowRange.preprocess(rM.rows());
    const auto processed_cols = rColRange.preprocess(rM.cols());
    return rM.block(processed_rows.start(), processed_cols.start(), processed_rows.size(), processed_cols.size());
}

template<class T, std::size_t N>
inline auto project(EigenBoundedVector<T, N>& rV, const boost::numeric::ublas::range& rRange)
{
    const auto processed = rRange.preprocess(rV.size());
    return rV.segment(processed.start(), processed.size());
}

template<class T, std::size_t N>
inline auto project(const EigenBoundedVector<T, N>& rV, const boost::numeric::ublas::range& rRange)
{
    const auto processed = rRange.preprocess(rV.size());
    return rV.segment(processed.start(), processed.size());
}

template<class T, std::size_t N1, std::size_t N2>
inline auto project(EigenBoundedMatrix<T, N1, N2>& rM, const boost::numeric::ublas::range& rRowRange, const boost::numeric::ublas::range& rColRange)
{
    const auto processed_rows = rRowRange.preprocess(rM.rows());
    const auto processed_cols = rColRange.preprocess(rM.cols());
    return rM.block(processed_rows.start(), processed_cols.start(), processed_rows.size(), processed_cols.size());
}

template<class T, std::size_t N1, std::size_t N2>
inline auto project(const EigenBoundedMatrix<T, N1, N2>& rM, const boost::numeric::ublas::range& rRowRange, const boost::numeric::ublas::range& rColRange)
{
    const auto processed_rows = rRowRange.preprocess(rM.rows());
    const auto processed_cols = rColRange.preprocess(rM.cols());
    return rM.block(processed_rows.start(), processed_cols.start(), processed_rows.size(), processed_cols.size());
}

///@}

} // namespace Kratos
