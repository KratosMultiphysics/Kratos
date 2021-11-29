//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_HELMHOLTZ_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_HELMHOLTZ_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
    KRATOS_DEFINE_VARIABLE( int, HELMHOLTZ_DIRECTION )
    KRATOS_DEFINE_VARIABLE( double, HELMHOLTZ_POISSON_RATIO )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_APPLICATION, HELMHOLTZ_VARS )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_APPLICATION, HELMHOLTZ_SOURCE )
}

#endif	/* KRATOS_HELMHOLTZ_APPLICATION_VARIABLES_H_INCLUDED */
