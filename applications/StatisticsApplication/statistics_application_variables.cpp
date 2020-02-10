//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#include "statistics_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_SUM)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_MEAN)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_VARIANCE)

KRATOS_CREATE_VARIABLE(double, VELOCITY_NORM)
KRATOS_CREATE_VARIABLE(double, PRESSURE_NORM)
KRATOS_CREATE_VARIABLE(double, PRESSURE_SUM)
KRATOS_CREATE_VARIABLE(double, PRESSURE_MEAN)
KRATOS_CREATE_VARIABLE(double, PRESSURE_VARIANCE)

} // namespace Kratos
