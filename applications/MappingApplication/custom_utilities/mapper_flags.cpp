//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

#if !defined(KRATOS_MAPPER_FLAGS_CPP_INCLUDED )
#define  KRATOS_MAPPER_FLAGS_CPP_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "mapper_flags.h"


namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, SWAP_SIGN,                0 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, ADD_VALUES,               1 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, POINT_WISE_VALUES,        2 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, CONSERVATIVE,             3 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, REMESHED,                 4 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, NON_CONFORMING_INTERFACE, 5 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, INVERSE_DIRECTION,        6 );
    KRATOS_CREATE_LOCAL_FLAG( MapperFlags, INTERPOLATE_VALUES,       7 );
}  // namespace Kratos.

#endif // KRATOS_MAPPER_FLAGS_CPP_INCLUDED  defined
