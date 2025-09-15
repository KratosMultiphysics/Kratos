//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#include "statistics_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_SUM)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_MEAN)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_VARIANCE)

KRATOS_CREATE_VARIABLE(double, VECTOR_3D_NORM)
KRATOS_CREATE_VARIABLE(double, SCALAR_NORM)
KRATOS_CREATE_VARIABLE(double, SCALAR_SUM)
KRATOS_CREATE_VARIABLE(double, SCALAR_MEAN)
KRATOS_CREATE_VARIABLE(double, SCALAR_VARIANCE)

} // namespace Kratos
