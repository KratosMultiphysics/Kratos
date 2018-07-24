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

#include "iga_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, NURBS_CONTROLPOINT_WEIGHT)

KRATOS_CREATE_VARIABLE(Vector, COORDINATES)
KRATOS_CREATE_VARIABLE(Vector, TANGENTS)

KRATOS_CREATE_VARIABLE(double, CROSS_AREA)
KRATOS_CREATE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_CREATE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)

KRATOS_CREATE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_BETA)

} // namespace Kratos
