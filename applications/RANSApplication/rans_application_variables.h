//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
    // incompressible potential flow specific variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, VELOCITY_POTENTIAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, PRESSURE_POTENTIAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_INLET )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_OUTLET )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_STRUCTURE )

} // namespace Kratos

#endif /* KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED */
