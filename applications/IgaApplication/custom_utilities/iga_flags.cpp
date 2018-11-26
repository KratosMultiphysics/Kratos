//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_IGA_FLAGS_CPP_INCLUDED )
#define  KRATOS_IGA_FLAGS_CPP_INCLUDED

// System includes

// External includes

// Project includes
#include "iga_flags.h"


namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_DISPLACEMENT_X,        0 );
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_DISPLACEMENT_Y,        1 );
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_DISPLACEMENT_Z,        2 );
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_ROTATION_X,            3 );
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_ROTATION_Y,            4 );
    KRATOS_CREATE_LOCAL_FLAG(IGAFlags, FIX_ROTATION_Z,            5 );
}  // namespace Kratos.

#endif // KRATOS_IGA_FLAGS_CPP_INCLUDED  defined