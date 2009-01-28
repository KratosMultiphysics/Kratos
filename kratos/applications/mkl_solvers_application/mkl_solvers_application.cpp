//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-07-23 14:55:54 $
//   Revision:            $Revision: 1.1 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "mkl_solvers_application.h"


namespace Kratos
{
	
	void KratosMKLSolversApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosMKLSolversApplication... " << std::endl;

	}


}  // namespace Kratos.


