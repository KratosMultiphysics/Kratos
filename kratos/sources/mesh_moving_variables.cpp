//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mesh_moving_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_ACCELERATION)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS)

    KRATOS_CREATE_VARIABLE(int, LAPLACIAN_DIRECTION)
    KRATOS_CREATE_VARIABLE(double, MESH_POISSON_RATIO)


void KratosApplication::RegisterALEVariables()
{
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_ACCELERATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS)

    KRATOS_REGISTER_VARIABLE(LAPLACIAN_DIRECTION)
    KRATOS_REGISTER_VARIABLE(MESH_POISSON_RATIO)

}

}  // namespace Kratos.
