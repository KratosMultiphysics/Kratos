//   
//   Project Name:        Kratos Wind Turbine Application
//   Last Modified by:    $Author:  $
//   Date:                $Date: 20120507 18:08:58 $
//   Revision:            $Revision: 0.1 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "wind_turbine_application.h"
#include "includes/variables.h"

namespace Kratos
{
        KRATOS_CREATE_VARIABLE(int, AUX_ID);
        KRATOS_CREATE_VARIABLE(int, AUX_INTERF_ID);
        KRATOS_CREATE_VARIABLE(int, AUX_BASE_ID);
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR);
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

	/// default constructor
 	KratosWindTurbineApplication::KratosWindTurbineApplication()
 	{}

 	void KratosWindTurbineApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosWindTurbineApplication... " << std::endl;
 
                KRATOS_REGISTER_VARIABLE( AUX_ID );
                KRATOS_REGISTER_VARIABLE( AUX_INTERF_ID );
                KRATOS_REGISTER_VARIABLE( AUX_BASE_ID );
// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR );
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 
 	}

}  // namespace Kratos.

