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
#include<omp.h>

// External includes


// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
template<class TDataType>
inline void AtomicAdd(TDataType& target, const TDataType& value ) {
    #pragma omp atomic
    target += value;
}

template<class TVectorType1, class TVectorType2>
inline void AtomicAdd(TVectorType1& target, const TVectorType2& value ) {
    for(unsigned int i=0: i<target.size(); ++i){
       AtomicAdd(target[i], value[i];)
    }
}

template<class TDataType>
inline void AtomicAdd(TDataType& target, const TDataType& value ) {
    #pragma omp atomic
    target += value;
}

template<class TVectorType1, class TVectorType2>
inline void AtomicSub(TVectorType1& target, const TVectorType2& value ) {
    for(unsigned int i=0: i<target.size(); ++i){
       AtomicSub(target[i], value[i];)
    }
}

template<class TDataType>
inline void AtomicAssign(TDataType& target, const TDataType& value) {
    #pragma omp atomic write
    target = value;
}

}  // namespace Kratos.

#endif // KRATOS_ATOMIC_UTILITIES_H_INCLUDED  defined
