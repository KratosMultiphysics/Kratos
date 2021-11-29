//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

#include "helmholtz_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE( int, HELMHOLTZ_DIRECTION )
    KRATOS_CREATE_VARIABLE( double, HELMHOLTZ_POISSON_RATIO )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_VARS )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_SOURCE )
}
