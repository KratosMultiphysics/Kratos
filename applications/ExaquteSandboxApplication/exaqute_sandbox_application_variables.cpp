//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

#include "exaqute_sandbox_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, VELOCITY_H1_SEMINORM )
KRATOS_CREATE_VARIABLE( double, AVERAGED_DIVERGENCE )
KRATOS_CREATE_VARIABLE( double, DIVERGENCE )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_NOISE )
KRATOS_CREATE_VARIABLE( double, AVERAGED_PRESSURE )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AVERAGED_VELOCITY )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DRAG_FORCE_CENTER_AUX)

}
