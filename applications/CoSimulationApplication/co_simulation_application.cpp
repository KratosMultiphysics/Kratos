//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "co_simulation_application.h"
#include "co_simulation_application_variables.h"


namespace Kratos {

KratosCoSimulationApplication::KratosCoSimulationApplication() {}

void KratosCoSimulationApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosCoSimulationApplication... " << std::endl;


}
}  // namespace Kratos.
