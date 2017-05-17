//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//

#include "pure_diffusion_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, POINT_HEAT_SOURCE)
KRATOS_CREATE_VARIABLE( double, HEIGHT)
KRATOS_CREATE_VARIABLE( double, PROJECTED_HEIGHT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)

}
