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


// System includes


// External includes


// Project includes
#include "droplet_dynamics_application.h"
#include "droplet_dynamics_application_variables.h"


namespace Kratos {

KratosDropletDynamicsApplication::KratosDropletDynamicsApplication():
    KratosApplication("DropletDynamicsApplication")
    {}

void KratosDropletDynamicsApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosDropletDynamicsApplication..." << std::endl;

}
}  // namespace Kratos.
