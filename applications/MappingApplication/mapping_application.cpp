//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
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

  std::cout << "    KRATOS ______  ___                      _____  "                          << std::endl;
  std::cout << "           ___   |/  /_____ ___________________(_)_____________ _  "          << std::endl;
  std::cout << "           __  /|_/ /_  __ `/__  __ \\__  __ \\_  /__  __ \\_  __ `/  "       << std::endl;
  std::cout << "           _  /  / / / /_/ /__  /_/ /_  /_/ /  / _  / / /  /_/ /  "           << std::endl;
  std::cout << "           /_/  /_/  \\__,_/ _  .___/_  .___//_/  /_/ /_/_\\__, /  "          << std::endl;
  std::cout << "                            /_/     /_/                 /____/ Application"   << std::endl;

 	std::cout << "Initializing KratosMappingApplication... " << std::endl;

  KRATOS_REGISTER_VARIABLE( NEIGHBOR_RANK )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEIGHBOR_COORDINATES )

}
}  // namespace Kratos.