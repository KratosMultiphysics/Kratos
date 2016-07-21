//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE(int,PATCH_INDEX)
KRATOS_CREATE_VARIABLE(double,TAUONE)
KRATOS_CREATE_VARIABLE(double,TAUTWO)
KRATOS_CREATE_VARIABLE(double,PRESSURE_MASSMATRIX_COEFFICIENT)

//KRATOS_CREATE_VARIABLE(double,Y_WALL)
KRATOS_CREATE_VARIABLE(double,SUBSCALE_PRESSURE)
KRATOS_CREATE_VARIABLE(double, C_DES)
//    KRATOS_CREATE_VARIABLE(double, C_SMAGORINSKY)


KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VORTICITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY)

// Non-Newtonian constitutive relations
KRATOS_CREATE_VARIABLE(double, REGULARIZATION_COEFFICIENT)

KRATOS_CREATE_VARIABLE(double, BINGHAM_SMOOTHER)
KRATOS_CREATE_VARIABLE(double, GEL_STRENGTH )

// Q-Criterion (for vortex visualization)
KRATOS_CREATE_VARIABLE(double, Q_VALUE)

// Vorticity
KRATOS_CREATE_VARIABLE(double, VORTICITY_MAGNITUDE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT)

// For swimming DEM
KRATOS_CREATE_VARIABLE(Vector, NODAL_WEIGHTS)
}
