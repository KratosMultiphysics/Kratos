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
    // Shallow water variables
    KRATOS_CREATE_VARIABLE( double, HEIGHT)                                 // Main variable
    KRATOS_CREATE_VARIABLE( double, BATHYMETRY)                             // Geometric definition of the problem
    KRATOS_CREATE_VARIABLE( double, RAIN)                                   // Source term
    KRATOS_CREATE_VARIABLE( double, FREE_SURFACE_ELEVATION)                 // Free surface elevation from z=0 (HEIGHT = FREE_SURFACE - BATHYMETRY)
    KRATOS_CREATE_VARIABLE( double, MANNING)                                // Friction coefficient

    // Specific variableS for PFEM2
    KRATOS_CREATE_VARIABLE( double, MEAN_SIZE)
    KRATOS_CREATE_VARIABLE( double, MEAN_VEL_OVER_ELEM_SIZE)
    KRATOS_CREATE_VARIABLE( double, PROJECTED_SCALAR1)
    KRATOS_CREATE_VARIABLE( double, DELTA_SCALAR1)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VECTOR1)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VECTOR1)

    // Units conversion
    KRATOS_CREATE_VARIABLE( double, TIME_UNIT_CONVERTER)
    KRATOS_CREATE_VARIABLE( double, WATER_HEIGHT_UNIT_CONVERTER)
}
