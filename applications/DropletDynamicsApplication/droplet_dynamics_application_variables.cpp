//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

#include "droplet_dynamics_application_variables.h"

namespace Kratos
{
    // External interfacial force, e.g. for including the electromagentic coupling
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXT_INT_FORCE)
}
