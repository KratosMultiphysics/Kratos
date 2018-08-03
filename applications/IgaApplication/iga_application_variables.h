/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
*/

#if !defined(KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED)
#define  KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{

KRATOS_DEFINE_VARIABLE(double, NURBS_CONTROL_POINT_WEIGHT)

KRATOS_DEFINE_VARIABLE(Vector, COORDINATES)
KRATOS_DEFINE_VARIABLE(Vector, TANGENTS)

KRATOS_DEFINE_VARIABLE(double, CROSS_AREA)
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)

KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_BETA)

} // namespace Kratos

#endif // !defined(KRATOS_IGA_APPLICATION_VARIABLES_H_INCLUDED)