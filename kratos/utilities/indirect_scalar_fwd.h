//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//  This header exists to avoid including indirect_scalar.h in variables.h.
//

#if !defined(KRATOS_INDIRECT_SCALAR_FWD_H_INCLUDED)
#define  KRATOS_INDIRECT_SCALAR_FWD_H_INCLUDED

#include <type_traits>

namespace Kratos
{

template <class T, typename S = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class IndirectScalar;

}  // namespace Kratos.

#endif // KRATOS_INDIRECT_SCALAR_FWD_H_INCLUDED  defined 
