//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.19 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "mpi_search_application.h"

namespace Kratos
{
  
	KratosMPISearchApplication::KratosMPISearchApplication()
	{}
 	
	void KratosMPISearchApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
		std::cout << "Initializing Kratos MPISearchApplication... " << std::endl;
	}
	
}  // namespace Kratos.


