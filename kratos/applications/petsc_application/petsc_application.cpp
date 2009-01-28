//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-04-21 13:13:16 $
//   Revision:            $Revision: 1.2 $
//
// 



// System includes


// External includes 

// Project includes
#include "includes/define.h"
#include "petsc_application.h"
#include "petsc_solver.h"


namespace Kratos
{
	
	void KratosPetscApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosPetscApplication... " << std::endl;
		int argc = 2;
		char** args;
		args = new char*[argc];
		for(int i = 0 ; i < argc ; i++)
		{
			args[i] = new char[255];
		}
		strcpy(args[0], "-ksp_type");
		strcpy(args[1], "bcgs");
//  		strcpy(args[2], "-ksp_max_it 100000");
// 		strcpy(args[2], "-ksp_monitor");

		char* help = 0;
		PetscInitialize(&argc,&args,(char *)0,help);

	}


}  // namespace Kratos.



  // workaround ubuntu bug	
// #define QUEUESTRINGSIZE 8192
// typedef struct _PrintfQueue *PrintfQueue;
// struct _PrintfQueue {
//   char        string[QUEUESTRINGSIZE];
//   PrintfQueue next;
// };
// PrintfQueue queue       = 0,queuebase = 0;
// int         queuelength = 0;
// FILE        *queuefile  = PETSC_NULL;
