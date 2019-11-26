// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "co_simulation_application.h"
#include "co_simulation_application_variables.h"

namespace Kratos {

KratosCoSimulationApplication::KratosCoSimulationApplication():
    KratosApplication("CoSimulationApplication")
    {}

void KratosCoSimulationApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
	KRATOS_INFO("") << "    KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ __\n"
	                << "           | |   / _ \\___ \\| | '_ ` _ \\| | | | |/ _` | __| |/ _ \\| '_ \\\n"
	                << "           | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |\n"
	                << "            \\____\\___/____/|_|_| |_| |_|\\__,_|_|\\__,_|\\__|_|\\___/|_| |_|\n"
                    << "Initializing KratosCoSimulationApplication..." << std::endl;


    KRATOS_REGISTER_VARIABLE(SCALAR_DISPLACEMENT);
    KRATOS_REGISTER_VARIABLE(SCALAR_ROOT_POINT_DISPLACEMENT);
    KRATOS_REGISTER_VARIABLE(SCALAR_REACTION);
    KRATOS_REGISTER_VARIABLE(SCALAR_FORCE);
    KRATOS_REGISTER_VARIABLE(SCALAR_VOLUME_ACCELERATION);


}
}  // namespace Kratos.
