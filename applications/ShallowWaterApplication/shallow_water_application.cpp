//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Miguel Masó Sotomayor $
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
#include "shallow_water_application.h"
#include "includes/variables.h"


namespace Kratos
{

 	KratosShallowWaterApplication::KratosShallowWaterApplication()
 	{}
 	
 	void KratosShallowWaterApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosShallowWaterApplication... " << std::endl; 

 	}

}  // namespace Kratos.


