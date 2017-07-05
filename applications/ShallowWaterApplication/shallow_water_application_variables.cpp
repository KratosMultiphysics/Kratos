//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//

#include "shallow_water_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, BATHYMETRY)

KRATOS_CREATE_VARIABLE( double, HEIGHT)
KRATOS_CREATE_VARIABLE( double, PROJECTED_HEIGHT)
KRATOS_CREATE_VARIABLE( double, DELTA_HEIGHT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_MOMENTUM)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_MOMENTUM)

KRATOS_CREATE_VARIABLE( double, MEAN_SIZE)
KRATOS_CREATE_VARIABLE( double, MEAN_VEL_OVER_ELEM_SIZE)

KRATOS_CREATE_VARIABLE( double, SCALAR_VELOCITY_X)
KRATOS_CREATE_VARIABLE( double, SCALAR_VELOCITY_Y)
KRATOS_CREATE_VARIABLE( double, SCALAR_PROJECTED_VELOCITY_X)
KRATOS_CREATE_VARIABLE( double, SCALAR_PROJECTED_VELOCITY_Y)
}
