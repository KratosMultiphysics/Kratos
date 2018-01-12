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
#include "stabilized_cfd_application.h"
#include "stabilized_cfd_application_variables.h"


namespace Kratos {

KratosStabilizedCFDApplication::KratosStabilizedCFDApplication():
    KratosApplication("StabilizedCFDApplication")
    {}

void KratosStabilizedCFDApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosStabilizedCFDApplication... " << std::endl;

  KRATOS_REGISTER_VARIABLE( FIC_BETA )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DIRECTIONAL_BETA )
  KRATOS_REGISTER_VARIABLE( RECORDED_STEPS )
  KRATOS_REGISTER_VARIABLE( MEAN_KINETIC_ENERGY )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MEAN_VELOCITY )
  KRATOS_REGISTER_VARIABLE( MEAN_PRESSURE )
  KRATOS_REGISTER_VARIABLE( VELOCITY_COVARIANCES )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_STATISTICS )
  KRATOS_REGISTER_VARIABLE( TRACE_XI )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DIV_XI )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION_RHS )
  KRATOS_REGISTER_VARIABLE( MASS_PROJECTION )
  KRATOS_REGISTER_VARIABLE( MASS_PROJECTION_RHS )

}
}  // namespace Kratos.
