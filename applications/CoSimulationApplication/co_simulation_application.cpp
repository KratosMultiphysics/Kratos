//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "co_simulation_application.h"

namespace Kratos {

KratosCoSimulationApplication::KratosCoSimulationApplication():
    KratosApplication("CoSimulationApplication")
    {}

void KratosCoSimulationApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosCoSimulationApplication..." << std::endl;
}
}  // namespace Kratos.
