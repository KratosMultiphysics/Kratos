//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-07-23 14:46:21 $
//   Revision:            $Revision: 1.2 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "externalsolvers_application.h"


namespace Kratos
{
	
	void KratosExternalSolversApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosExternalSolversApplication... " << std::endl;

	}


}  // namespace Kratos.


