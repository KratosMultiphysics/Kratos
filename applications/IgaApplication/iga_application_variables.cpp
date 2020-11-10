/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#include "iga_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, CROSS_AREA)
KRATOS_CREATE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_CREATE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_BETA)

//Load Condition Variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DEAD_LOAD)

//Penalty Variables
KRATOS_CREATE_VARIABLE(double, PENALTY_FACTOR)

} // namespace Kratos
