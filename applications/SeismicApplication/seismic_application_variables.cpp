//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mahmoud Zidan
//

#include "seismic_application_variables.h"

namespace Kratos
{

// Fiber Beam-Column Element Variables
KRATOS_CREATE_VARIABLE(double, ELEMENT_LOOP_TOLERANCE)
KRATOS_CREATE_VARIABLE(int, MAX_EQUILIBRIUM_ITERATIONS)
KRATOS_CREATE_VARIABLE(int, NUMBER_OF_SECTIONS)
KRATOS_CREATE_VARIABLE(Matrix, CONCRETE_FIBERS_DATA)
KRATOS_CREATE_VARIABLE(Matrix, STEEL_FIBERS_DATA)
KRATOS_CREATE_VARIABLE(double, BEAM_WIDTH)
KRATOS_CREATE_VARIABLE(double, BEAM_HEIGHT)

// Fiber Beam-Column Constitutive parameters
KRATOS_CREATE_VARIABLE(double, CONCRETE_YIELD_STRENGTH)
KRATOS_CREATE_VARIABLE(double, CONCRETE_YIELD_STRAIN)
KRATOS_CREATE_VARIABLE(double, CONCRETE_CRUSHING_STRAIN)
KRATOS_CREATE_VARIABLE(double, STEEL_YOUNGS_MODULUS)
KRATOS_CREATE_VARIABLE(double, STEEL_HARDENING_RATIO)
KRATOS_CREATE_VARIABLE(double, STEEL_TRANSITION_VARIABLE)
KRATOS_CREATE_VARIABLE(double, STEEL_YIELD_STRENGTH)
KRATOS_CREATE_VARIABLE(double, STEEL_A1_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, STEEL_A2_COEFFICIENT)

}
