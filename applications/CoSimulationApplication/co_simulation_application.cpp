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
	std::cout  << "KRATOS / ___|  ___   ___ (_) _ __   _   _| | __ _| |_(_)            " << std::endl;
	std::cout  << "       | |    / _ \\/ __| | '_ ` _ \\| | | | |/ _` | __| |/ _ \\| '_ \\ " << std::endl;
	std::cout  << "       | |__ | (_) |\\__ \\| | | | | | |_| | | (_| | |_| | (_) | | | |" << std::endl;
	std::cout  << "       \\____| \\___/___/_/| | | | |\\__,_|_|__\\__,_|\\__|_|\\___/|_| |_| Application" << std::endl;
    std::cout << "Initializing KratosCoSimulationApplication... " << std::endl;
}
} // namespace Kratos.
