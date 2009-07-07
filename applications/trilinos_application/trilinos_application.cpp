//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-11-10 14:24:45 $
//   Revision:            $Revision: 1.1 $
//
// 



// System includes


// External includes 

// Project includes
#include "includes/define.h"
#include "trilinos_application.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_DataAccess.h"

// #include "Epetra_MpiComm.h"

//#include "trilinos_solver.h"

namespace Kratos
{
		
	KRATOS_CREATE_VARIABLE( bool, IS_INACTIVE )


	void KratosTrilinosApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosTrilinosApplication... " << std::endl;
		KRATOS_REGISTER_VARIABLE( IS_INACTIVE )
	}


}  // namespace Kratos.


