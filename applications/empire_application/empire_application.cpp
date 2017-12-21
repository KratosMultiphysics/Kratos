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
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "empire_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosEmpireApplication::KratosEmpireApplication() : KratosApplication("EmpireApplication") {}
 	
 	void KratosEmpireApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();

    std::stringstream banner;
    banner << "    KRATOS  _____ __  __ ____ ___ ____  _____ "               << std::endl;
    banner << "           | ____|  \\/  |  _ \\_ _|  _ \\| ____|"            << std::endl;
    banner << "           |  _| | |\\/| | |_) | || |_) |  _| "               << std::endl;
    banner << "           | |___| |  | |  __/| ||  _ <| |___ "               << std::endl;
    banner << "           |_____|_|  |_|_|  |___|_| \\_\\_____| Application" << std::endl;

    banner << "Initializing KratosEmpireApplication... " << std::endl;
    
    std::cout << banner.str();
 
// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 	}

}  // namespace Kratos.
                                    
