//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//

#include "constitutive_laws_application_variables.h"

namespace Kratos
{
// Constitutive laws variables
KRATOS_CREATE_VARIABLE(Vector, HIGH_CYCLE_FATIGUE_COEFFICIENTS)
KRATOS_CREATE_VARIABLE(double, FATIGUE_REDUCTION_FACTOR)
KRATOS_CREATE_VARIABLE(int, NUMBER_OF_CYCLES)
KRATOS_CREATE_VARIABLE(int, LOCAL_NUMBER_OF_CYCLES)
KRATOS_CREATE_VARIABLE(double, WOHLER_STRESS)
KRATOS_CREATE_VARIABLE(double, REVERSION_FACTOR_RELATIVE_ERROR)
KRATOS_CREATE_VARIABLE(double, MAX_STRESS_RELATIVE_ERROR)
KRATOS_CREATE_VARIABLE(double, MAX_STRESS)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_STRESS)
KRATOS_CREATE_VARIABLE(bool, CYCLE_INDICATOR)
KRATOS_CREATE_VARIABLE(double, CYCLES_TO_FAILURE)
KRATOS_CREATE_VARIABLE(double, TIME_INCREMENT)
KRATOS_CREATE_VARIABLE(bool, DAMAGE_ACTIVATION)
KRATOS_CREATE_VARIABLE(double, PREVIOUS_CYCLE);
KRATOS_CREATE_VARIABLE(double, CYCLE_PERIOD)
KRATOS_CREATE_VARIABLE(bool, ADVANCE_STRATEGY_APPLIED);

}
