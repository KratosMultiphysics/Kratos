//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Denis Demidov
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <atomic>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

#if !defined(__cpp_lib_atomic_ref) && defined(KRATOS_SMP_CXX11)
#include <boost/atomic/atomic_ref.hpp>
#endif

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

namespace Kratos {

#if defined(KRATOS_SMP_CXX11) || defined(KRATOS_COMPILED_IN_OS)
    #if defined(__cpp_lib_atomic_ref) // C++20
        template <class T>
        using AtomicRef = std::atomic_ref<T>;
    #else
        template <class T>
        using AtomicRef = boost::atomic_ref<T>;
    #endif //__cpp_lib_atomic_ref
#endif // KRATOS_SMP_CXX11

///@addtogroup KratosCore
/**
 * collection of utilities for atomic updates of simple types. (essentially mimics the omp atomic)
 */

/**
 * @param target variable being atomically updated by doing target += value
 * @param value value being added
 */
template<class TDataType>
inline void AtomicAdd(TDataType& target, const TDataType& value)
{
#if defined(KRATOS_SMP_CXX11) || defined(KRATOS_COMPILED_IN_OS)
    AtomicRef<TDataType>{target} += value;
#elif defined(KRATOS_SMP_CXX11)
    #pragma omp atomic
    target += value;
#else
    target += value;
#endif
}

/**
 * @param target variable being atomically updated by doing target += value
 * @param value value being added
 * Specialization for array_1d
 * Note that the update is not really atomic, but rather is done component by component
 */
template <class TDataType, std::size_t ArraySize>
inline void AtomicAdd(array_1d<TDataType,ArraySize>& target, const array_1d<TDataType,ArraySize>& value)
{
    for (std::size_t i=0; i<ArraySize; ++i) {
        AtomicAdd(target[i], value[i]);
    }
}

/**
 * @param target vector variable being atomically updated by doing target += value
 * @param value vector value being added
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicAddVector(TVectorType1& target, const TVectorType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicAddVector- Sizes are: " << target.size() << " for target and " << value.size() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size(); ++i) {
       AtomicAdd(target[i], value[i]);
    }
}

/**
 * @param target matrix variable being atomically updated by doing target -= value
 * @param value matrix value being subtracted
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TMatrixType1, class TMatrixType2>
inline void AtomicAddMatrix(TMatrixType1& target, const TMatrixType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size1() != value.size1() || target.size2() != value.size2()) << "matrix size mismatch in matrix AtomicAddMatrix- Sizes are: " << target.size1() << "x" << target.size2() << " for target and " << value.size1() << "x" << value.size2() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size1(); ++i) {
        for(std::size_t j=0; j<target.size2(); ++j) {
            AtomicAdd(target(i,j), value(i,j));
        }
    }
}

/**
 * @param target vector variable being atomically updated by doing target -= value
 * @param value vector value being subtracted
 */
template<class TDataType>
inline void AtomicSub(TDataType& target, const TDataType& value)
{
    AtomicAdd(target, -value);
}

/**
 * @param target variable being atomically updated by doing target -= value
 * @param value value being subtracted
 * Specialization for array_1d
 * Note that the update is not really atomic, but rather is done component by component
 */
template <class TDataType, std::size_t ArraySize>
inline void AtomicSub(array_1d<TDataType,ArraySize>& target, const array_1d<TDataType,ArraySize>& value)
{
    for(std::size_t i=0; i<ArraySize; ++i) {
        AtomicSub(target[i], value[i]);
    }
}

/**
 * @param target vector variable being atomically updated by doing target -= value
 * @param value vector value being subtracted
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicSubVector(TVectorType1& target, const TVectorType2& value) {
    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicSubVector- Sizes are: " << target.size() << " for target and " << value.size() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size(); ++i) {
       AtomicSub(target[i], value[i]);
    }
}

/**
 * @param target matrix variable being atomically updated by doing target -= value
 * @param value matrix value being subtracted
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TMatrixType1, class TMatrixType2>
inline void AtomicSubMatrix(TMatrixType1& target, const TMatrixType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size1() != value.size1() || target.size2() != value.size2()) << "matrix size mismatch in matrix AtomicSubMatrix- Sizes are: " << target.size1() << "x" << target.size2() << " for target and " << value.size1() << "x" << value.size2() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size1(); ++i) {
        for(std::size_t j=0; j<target.size2(); ++j) {
            AtomicSub(target(i,j), value(i,j));
        }
    }
}

/** @param ref scalar variable being atomically updated by doing ref *= factor
 * @param factor vector value being multiplied
 */
template<class T>
inline T atomic_mul(std::atomic_ref<T>& ref, T factor)
{
    T old = ref.load(std::memory_order_relaxed);
    T result;
    do {
        result = old * result;
    } while (!ref.compare_exchange_weak(old, result, std::memory_order_acq_rel, std::memory_order_relaxed));

    return result;
}

/** @param target vector variable being atomically updated by doing target *= value
 * @param value vector value being multiplied
 * Note that "=*" compound operator is --> **NOT** <-- defined for multiplication in boost::atomic//std::atomic,
 * so we need to implement it ourselves
 */
template<class TDataType>
inline void AtomicMult(TDataType& target, const TDataType& value)
{
#if defined(KRATOS_SMP_CXX11) || defined(KRATOS_COMPILED_IN_OS)
    AtomicRef<TDataType> at_ref{target};
    atomic_mul(at_ref, value);
#elif defined(KRATOS_SMP_OPENMP)
    #pragma omp atomic
    target *= value;
#else
    target *= value;
#endif
}

/** @param target variable being atomically updated by doing target *= value
 * @param value value being multiplied
 * Specialization for array_1d
 * Note that the update is not really atomic, but rather is done component by component
 */
template <class TDataType, std::size_t ArraySize>
inline void AtomicMult(array_1d<TDataType,ArraySize>& target, const array_1d<TDataType,ArraySize>& value)
{
    for(std::size_t i=0; i<ArraySize; ++i) {
        AtomicMult(target[i], value[i]);
    }
}

/** @param target vector variable being atomically updated by doing target *= value
 * @param value vector value being multiplied
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicMultVector(TVectorType1& target, const TVectorType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicMultVector- Sizes are: " << target.size() << " for target and " << value.size() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size(); ++i) {
       AtomicMult(target[i], value[i]);
    }
}

/**
 * @param target matrix variable being atomically updated by doing target *= value
 * @param value matrix value being multiplied
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TMatrixType1, class TMatrixType2>
inline void AtomicMultMatrix(TMatrixType1& target, const TMatrixType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size1() != value.size1() || target.size2() != value.size2()) << "matrix size mismatch in matrix AtomicMultMatrix- Sizes are: " << target.size1() << "x" << target.size2() << " for target and " << value.size1() << "x" << value.size2() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size1(); ++i) {
        for(std::size_t j=0; j<target.size2(); ++j) {
            AtomicMult(target(i,j), value(i,j));
        }
    }
}

/** @param target vector variable being atomically updated by doing target *= 1.0/value
 * @param value vector value being divided
 */
template<class TDataType>
inline void AtomicDiv(TDataType& target, const TDataType& value)
{
    AtomicMult(target, 1.0/value);
}

/** @param target variable being atomically updated by doing target *= 1.0/value
 * @param value value being divided
 * Specialization for array_1d
 * Note that the update is not really atomic, but rather is done component by component
 */
template <class TDataType, std::size_t ArraySize>
inline void AtomicDiv(array_1d<TDataType,ArraySize>& target, const array_1d<TDataType,ArraySize>& value)
{
    for(std::size_t i=0; i<ArraySize; ++i) {
        AtomicDiv(target[i], value[i]);
    }
}

/** @param target vector variable being atomically updated by doing target *= 1.0/value
 * @param value vector value being divided
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicDivVector(TVectorType1& target, const TVectorType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicDivVector- Sizes are: " << target.size() << " for target and " << value.size() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size(); ++i) {
       AtomicDiv(target[i], value[i]);
    }
}

/**
 * @param target matrix variable being atomically updated by doing target *= 1.0/value
 * @param value matrix value being divided
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TMatrixType1, class TMatrixType2>
inline void AtomicDivMatrix(TMatrixType1& target, const TMatrixType2& value)
{
    KRATOS_DEBUG_ERROR_IF(target.size1() != value.size1() || target.size2() != value.size2()) << "matrix size mismatch in matrix AtomicDivMatrix- Sizes are: " << target.size1() << "x" << target.size2() << " for target and " << value.size1() << "x" << value.size2() << " for value " << std::endl;

    for(std::size_t i=0; i<target.size1(); ++i) {
        for(std::size_t j=0; j<target.size2(); ++j) {
            AtomicDiv(target(i,j), value(i,j));
        }
    }
}

}  // namespace Kratos.
