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
#include "mapping_application.h"
#include "mapping_application_variables.h"


namespace Kratos {

KratosMappingApplication::KratosMappingApplication() {}

void KratosMappingApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosMappingApplication... " << std::endl;


}
}  // namespace Kratos.
