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
KRATOS_CREATE_VARIABLE( double, BATHYMETRY)                             // Geometric definition of the problem
KRATOS_CREATE_VARIABLE( double, RAIN)                                   // Source term

KRATOS_CREATE_VARIABLE( double, HEIGHT)                                 // Main variable
KRATOS_CREATE_VARIABLE( double, PROJECTED_HEIGHT)                       // Convected variable
KRATOS_CREATE_VARIABLE( double, DELTA_HEIGHT)                           // Variable to update particles
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_MOMENTUM)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_MOMENTUM)

KRATOS_CREATE_VARIABLE( double, MEAN_SIZE)                              // Specific variable for PFEM2
KRATOS_CREATE_VARIABLE( double, MEAN_VEL_OVER_ELEM_SIZE)                // Specific variable for PFEM2

KRATOS_CREATE_VARIABLE( double, TIME_UNIT_CONVERTER)
KRATOS_CREATE_VARIABLE( double, WATER_HEIGHT_UNIT_CONVERTER)
KRATOS_CREATE_VARIABLE( double, HORIZONTAL_SCALE_UNIT_CONVERTER)
KRATOS_CREATE_VARIABLE( double, VERTICAL_SCALE_UNIT_CONVERTER)
}
