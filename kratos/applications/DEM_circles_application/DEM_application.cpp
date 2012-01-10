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


#include "includes/variables.h"
#include "DEM_application.h"


namespace Kratos
{
	/*
	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);*/


	KratosDEMApplication::KratosDEMApplication()
	{}
	
	void KratosDEMApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosDEMApplication... " << std::endl;
	}

}  // namespace Kratos.


