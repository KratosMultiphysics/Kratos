//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#include "constitutive_models_application_variables.h"

namespace Kratos
{

  //specific constitutive models variables must be CREATED here
  KRATOS_CREATE_VARIABLE(std::string, TEMPERATURE_VARIABLE)
  KRATOS_CREATE_VARIABLE(std::string, PRESSURE_VARIABLE)

  KRATOS_CREATE_VARIABLE(PropertiesLayout, PROPERTIES_LAYOUT)

  
  KRATOS_CREATE_VARIABLE( double, ALPHA )
   KRATOS_CREATE_VARIABLE( double, BETA )   
   KRATOS_CREATE_VARIABLE( double, MF )   
   KRATOS_CREATE_VARIABLE( double, CC )   
   KRATOS_CREATE_VARIABLE( double, MM )   
   KRATOS_CREATE_VARIABLE( double, RHOS )   
   KRATOS_CREATE_VARIABLE( double, RHOT )   
   KRATOS_CREATE_VARIABLE( double, KSIS )   
   KRATOS_CREATE_VARIABLE( double, RHOM )   
   KRATOS_CREATE_VARIABLE( double, KSIM )   
   KRATOS_CREATE_VARIABLE( double, PC0 )   
   KRATOS_CREATE_VARIABLE( double, CHIS )   
   KRATOS_CREATE_VARIABLE( double, CHIT )   

   KRATOS_CREATE_VARIABLE( double, VOID_RATIO )   
   KRATOS_CREATE_VARIABLE( double, PS )   
   KRATOS_CREATE_VARIABLE( double, PT )   
   KRATOS_CREATE_VARIABLE( double, PM )   
   KRATOS_CREATE_VARIABLE( double, PLASTIC_MULTIPLIER )   
}
