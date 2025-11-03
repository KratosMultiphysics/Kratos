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

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RESULTANT_FORCE);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RESULTANT_MOMENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(GLOBAL_DISPLACEMENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(GLOBAL_ROTATION);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(IMPOSED_FORCE);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(IMPOSED_MOMENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(IMPOSED_DISPLACEMENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(IMPOSED_ROTATION);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EFFECTIVE_FORCE);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EFFECTIVE_MOMENT);

    KRATOS_REGISTER_VARIABLE(NODES_ID_INDEX_MAP);
    KRATOS_REGISTER_VARIABLE(ELEMENTS_ID_INDEX_MAP);

    KRATOS_REGISTER_VARIABLE(COUPLING_ITERATION_NUMBER)

    KRATOS_REGISTER_VARIABLE(INTERFACE_EQUATION_ID)
    KRATOS_REGISTER_VARIABLE(EXPLICIT_EQUATION_ID)

    KRATOS_REGISTER_VARIABLE(MIDDLE_VELOCITY)

}
}  // namespace Kratos.
