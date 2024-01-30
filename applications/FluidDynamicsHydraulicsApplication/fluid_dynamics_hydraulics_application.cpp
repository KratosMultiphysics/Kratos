//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


// System includes


// External includes


// Project includes
#include "fluid_dynamics_hydraulics_application.h"
#include "fluid_dynamics_hydraulics_application_variables.h"


namespace Kratos {

KratosFluidDynamicsHydraulicsApplication::KratosFluidDynamicsHydraulicsApplication():
    KratosApplication("FluidDynamicsHydraulicsApplication")
    {}

void KratosFluidDynamicsHydraulicsApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosFluidDynamicsHydraulicsApplication..." << std::endl;


}

}  // namespace Kratos.
