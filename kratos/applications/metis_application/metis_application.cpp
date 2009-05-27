//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-07-09 15:09:17 $
//   Revision:            $Revision: 1.2 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"

#include "metis_application.h"
#include "includes/variables.h"
//#include "geometries/triangle_2d.h"


namespace Kratos
{

/*	KRATOS_CREATE_VARIABLE(double, NODAL_AREA)
	KRATOS_CREATE_VARIABLE(double, NODAL_H)
	KRATOS_CREATE_VARIABLE(double, IS_STRUCTURE)
	KRATOS_CREATE_VARIABLE(double, IS_FLUID)
	KRATOS_CREATE_VARIABLE(double, IS_BOUNDARY)
	KRATOS_CREATE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_CREATE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
*/
//	KRATOS_CREATE_VARIABLE(double, PRESSURE_OLD_IT)
//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX)

	KratosMetisApplication::KratosMetisApplication()	
	{}


	void KratosMetisApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing Kratos MetisApplication... " << std::endl;

	}


	
}  // namespace Kratos.


