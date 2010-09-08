//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "kratos_multibody_application.h"
#include "includes/variables.h"


namespace Kratos
{

  // Creating variables here
   KRATOS_CREATE_VARIABLE(double, IS_BODY )

      void KratosMultibodyApplication::Register()
	{
	  std::cout << "Initializig KratosMultibodyApplication... " << std::endl;

	  // Registering variables here
 	  KRATOS_REGISTER_VARIABLE(IS_BODY);

	// calling base class register to register Kratos components
	KratosApplication::Register();
	
	  // Registering elements and conditions here
// 	  KRATOS_REGISTER_ELEMENT("ApplicationElementName", msApplicationElement);
// 	  KRATOS_REGISTER_CONDITION("ApplicationConditionName", msApplicationCondition);

	}
      
      
  // Initializing static members
//   const ApplicationElement  KratosMultibodyApplication::msApplicationElement(arguments...);
//   const ApplicationCondition  KratosMultibodyApplication::msApplicationElement(arguments...); 

  
  
}  // namespace Kratos.


