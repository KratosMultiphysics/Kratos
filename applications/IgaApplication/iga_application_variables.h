/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED)
#define  KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

KRATOS_DEFINE_VARIABLE(double, NURBS_CONTROL_POINT_WEIGHT)

KRATOS_DEFINE_VARIABLE(Vector, COORDINATES)
KRATOS_DEFINE_VARIABLE(Vector, TANGENTS)

KRATOS_DEFINE_VARIABLE(double, CROSS_AREA)
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)
KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)

KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_BETA)

//Condition Variables
KRATOS_DEFINE_VARIABLE(Vector, POINT_LOAD)

KRATOS_DEFINE_VARIABLE(double, PENALTY_FACTOR)

} // namespace Kratos

#endif // !defined(KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED)