//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso
//

// System includes

// External includes

// Project includes
#include "meshing_flags.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG( MeshingFlags, REFINED,    0 );
KRATOS_CREATE_LOCAL_FLAG( MeshingFlags, TO_COARSE , 1 );
}  // namespace Kratos