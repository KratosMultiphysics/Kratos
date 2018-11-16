//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: Kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#if !defined(KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
    // Shallow water variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, HEIGHT)                               // Main variable
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, BATHYMETRY)                           // Geometric definition of the problem
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, RAIN)                                 // Source term
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, FREE_SURFACE_ELEVATION)               // Free surface elevation from z=0 (HEIGHT = FREE_SURFACE - BATHYMETRY)
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, MANNING)                              // Friction coefficient
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, TOPOGRAPHY_SLOPE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, TOPOGRAPHY_GRADIENT)

    // Specific variables for PFEM2
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, MEAN_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, MEAN_VEL_OVER_ELEM_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, PROJECTED_SCALAR1)
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, DELTA_SCALAR1)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, PROJECTED_VECTOR1)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, DELTA_VECTOR1)

    // Units conversion
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, TIME_UNIT_CONVERTER)
    KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, WATER_HEIGHT_UNIT_CONVERTER)
}

#endif	/* KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED */
