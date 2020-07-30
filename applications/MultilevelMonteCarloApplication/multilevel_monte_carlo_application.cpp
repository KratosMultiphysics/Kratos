//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes


// External includes


// Project includes
#include "multilevel_monte_carlo_application.h"
#include "multilevel_monte_carlo_application_variables.h"


namespace Kratos {

KratosMultilevelMonteCarloApplication::KratosMultilevelMonteCarloApplication():
    KratosApplication("MultilevelMonteCarloApplication")
    {}

void KratosMultilevelMonteCarloApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
	KRATOS_INFO("") <<
    " KRATOS __ __ _   __ __  __ " << std::endl <<
    "       |  V  | | |  V  |/ _|" << std::endl <<
    "       | \\_/ | |_| \\_/ | (_ " << std::endl <<
    "       |_| |_|___|_| |_|\\__|APPLICATION" << std::endl;

    KRATOS_REGISTER_VARIABLE( POWER_SUM_1 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_2 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_3 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_4 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_5 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_6 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_7 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_8 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_9 )
    KRATOS_REGISTER_VARIABLE( POWER_SUM_10 )
}
}  // namespace Kratos.
