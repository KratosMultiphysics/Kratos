//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "mapper_flags.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, SWAP_SIGN,                0 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, ADD_VALUES,               1 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, REMESHED,                 2 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, USE_TRANSPOSE,            3 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, ORIGIN_ONLY,              4 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, DESTINATION_ONLY,         5 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, TO_NON_HISTORICAL,        6 );
KRATOS_CREATE_LOCAL_FLAG( MapperFlags, FROM_NON_HISTORICAL,      7 );
}  // namespace Kratos.
