//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//


// System includes


// External includes


// Project includes
#include "contact_mechanics_application.h"
#include "contact_mechanics_application_variables.h"


namespace Kratos {

KratosContactMechanicsApplication::KratosContactMechanicsApplication() {}

void KratosContactMechanicsApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();

	std::cout << "             ___         _           _          " << std::endl;  
	std::cout << "     KRATOS / __|___ _ _| |_ __ _ __| |_          " << std::endl;
	std::cout << "           | (__/ _ \\ ' \\  _/ _` / _|  _|         " << std::endl;
	std::cout << "            \\___\\___/_||_\\__\\__,_\\__|\\__|MECHANICS" << std::endl;                     
 	std::cout << "Initializing KratosContactMechanicsApplication... " << std::endl;

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( STEP_ROTATION )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( STEP_DISPLACEMENT )

  KRATOS_REGISTER_VARIABLE( RESIDUAL_VECTOR )

  KRATOS_REGISTER_VARIABLE( NUMBER_OF_ACTIVE_CONTACTS )
  KRATOS_REGISTER_VARIABLE( NUMBER_OF_STICK_CONTACTS )
  KRATOS_REGISTER_VARIABLE( NUMBER_OF_SLIP_CONTACTS )
  KRATOS_REGISTER_VARIABLE( FRICTION_ACTIVE )
  KRATOS_REGISTER_VARIABLE( PENALTY_PARAMETER )
  KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL )
  KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL_REACTION )
  KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL )
  KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )
  KRATOS_REGISTER_VARIABLE( TAU_STAB )
  KRATOS_REGISTER_VARIABLE( MU_STATIC )
  KRATOS_REGISTER_VARIABLE( MU_DYNAMIC )

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )

  KRATOS_REGISTER_VARIABLE( CONTACT_ADHESION )
  KRATOS_REGISTER_VARIABLE( CONTACT_FRICTION_ANGLE )
  KRATOS_REGISTER_VARIABLE( TANGENTIAL_PENALTY_RATIO )
  KRATOS_REGISTER_VARIABLE( CONTACT_PLASTIC_SLIP )

}
}  // namespace Kratos.
