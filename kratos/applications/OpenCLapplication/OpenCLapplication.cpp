//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "OpenCLapplication.h"
#include "includes/variables.h"


namespace Kratos
{

	KratosOpenCLApplication::KratosOpenCLApplication()
	{}
 	
	void KratosOpenCLApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
		std::cout << "Initializing KratosTestApplication... " << std::endl;
	}

}  // namespace Kratos.


