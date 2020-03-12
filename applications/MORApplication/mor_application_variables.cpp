//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#include "mor_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, FREQUENCY )
KRATOS_CREATE_VARIABLE( double, ACOUSTIC_PRESSURE )
KRATOS_CREATE_VARIABLE( double, PRESSURE_GRADIENT )
KRATOS_CREATE_VARIABLE( double, ACOUSTIC_VELOCITY )
KRATOS_CREATE_VARIABLE( double, ACOUSTIC_PRESSURE_RESIDUAL)
KRATOS_CREATE_VARIABLE( double, ACOUSTIC_DISPLACEMENT )
}
