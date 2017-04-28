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

// This define must be HERE
#define DKRATOS_EXPORT_INTERFACE_2 1

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ale_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS)

void KratosApplication::RegisterALEVariables()
{
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS)
}

}  // namespace Kratos.

// This define must be HERE
#undef DKRATOS_EXPORT_INTERFACE_2
