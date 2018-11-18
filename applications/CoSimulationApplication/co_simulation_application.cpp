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
//

// System includes

// External includes

// Project includes
#include "co_simulation_application.h"

namespace Kratos
{

KratosCoSimulationApplication::KratosCoSimulationApplication() : KratosApplication("CoSimulationApplication") {}

void KratosCoSimulationApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
	KRATOS_INFO("") << "   KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ __                     " << std::endl
	                << "          | |   / _ \\___ \\| | '_ ` _ \\| | | | |/ _` | __| |/ _ \\| '_ \\               " << std::endl
	                << "          | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |                   " << std::endl
	                << "           \\____\\___/____/|_|_| |_| |_|\\__,_|_|\\__,_|\\__|_|\\___/|_| |_| Application " << std::endl
                    << "Initializing KratosCoSimulationApplication...										      " << std::endl;
}
} // namespace Kratos.
