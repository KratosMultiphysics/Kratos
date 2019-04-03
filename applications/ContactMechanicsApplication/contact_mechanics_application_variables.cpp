//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#include "contact_mechanics_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef FrictionLaw::Pointer  FrictionLawPointerType;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Create Variables

  KRATOS_CREATE_VARIABLE( std::string, FRICTION_LAW_NAME )
  KRATOS_CREATE_VARIABLE( FrictionLawPointerType, FRICTION_LAW )

  KRATOS_CREATE_VARIABLE( bool, FRICTION_ACTIVE )
  KRATOS_CREATE_VARIABLE( bool, HYDRAULIC )

  KRATOS_CREATE_VARIABLE( double, PENALTY_PARAMETER )
  KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_NORMAL )
  KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_NORMAL_REACTION )
  KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_TANGENTIAL )
  KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )

  KRATOS_CREATE_VARIABLE( double, TAU_STAB )
  KRATOS_CREATE_VARIABLE( double, MU_STATIC )
  KRATOS_CREATE_VARIABLE( double, MU_DYNAMIC )

  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WATER_CONTACT_FORCE )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )

  KRATOS_CREATE_VARIABLE( double, CONTACT_ADHESION )
  KRATOS_CREATE_VARIABLE( double, CONTACT_FRICTION_ANGLE )
  KRATOS_CREATE_VARIABLE( double, TANGENTIAL_PENALTY_RATIO )
  KRATOS_CREATE_VARIABLE( double, CONTACT_PLASTIC_SLIP )

  //thermal properties
  KRATOS_CREATE_VARIABLE( double, HEAT_CONDUCTIVITY )

  //solution
  KRATOS_CREATE_VARIABLE( int, SEGREGATED_STEP )
  KRATOS_CREATE_VARIABLE( double, CONTACT_STEP_TIME )

  ///@}
}
