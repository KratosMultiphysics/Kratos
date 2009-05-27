//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-09-28 12:56:44 $
//   Revision:            $Revision: 1.4 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "fsi_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
	KRATOS_CREATE_VARIABLE(double, AUX)
	KRATOS_CREATE_VARIABLE(double, FICTITIOUS_FLUID_DENSITY)
//	KRATOS_CREATE_VARIABLE(double, PRESSURE_OLD_IT)

//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RELAXED_DISP)
	
	void KratosFSIApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosFSIApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(AUX)
		KRATOS_REGISTER_VARIABLE(FICTITIOUS_FLUID_DENSITY)
//		KRATOS_REGISTER_VARIABLE(PRESSURE_OLD_IT)
//		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
//		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RELAXED_DISP)
	

	}

}  // namespace Kratos.


