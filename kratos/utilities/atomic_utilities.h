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
//

#if !defined(KRATOS_ATOMIC_UTILITIES_H_INCLUDED )
#define  KRATOS_ATOMIC_UTILITIES_H_INCLUDED


// System includes

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
/**
 * collection of utilities for atomic updates of simple types. (essentially mimics the omp atomic)
 */

/** @param target variable being atomically updated by doing target += value
 * @param value value being added
 */
template<class TDataType>
inline void AtomicAdd(TDataType& target, const TDataType& value ) {
    #pragma omp atomic
    target += value;
}

/** @param target vector variable being atomically updated by doing target += value
 * @param value vector value being added
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicAdd(TVectorType1& target, const TVectorType2& value ) {

    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicAdd- Sizes are: "
        << target.size() << " for target and " << value.size() << " for value " <<std::endl;
    for(unsigned int i=0; i<target.size(); ++i){
       AtomicAdd(target[i], value[i]);
    }
}

/** @param target vector variable being atomically updated by doing target -= value
 * @param value vector value being subtracted
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TDataType>
inline void AtomicSub(TDataType& target, const TDataType& value ) {
    #pragma omp atomic
    target -= value;
}

/** @param target vector variable being atomically updated by doing target -= value
 * @param value vector value being subtracted
 * Note that the update is not really atomic, but rather is done component by component
 */
template<class TVectorType1, class TVectorType2>
inline void AtomicSub(TVectorType1& target, const TVectorType2& value ) {
    KRATOS_DEBUG_ERROR_IF(target.size() != value.size()) << "vector size mismatch in vector AtomicSub- Sizes are: "
        << target.size() << " for target and " << value.size() << " for value " <<std::endl;
    for(unsigned int i=0; i<target.size(); ++i){
       AtomicSub(target[i], value[i]);
    }
}

}  // namespace Kratos.

#endif // KRATOS_ATOMIC_UTILITIES_H_INCLUDED  defined
