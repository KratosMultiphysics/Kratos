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
#include "swimming_DEM_application.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "../DEM_application/DEM_application.h"

namespace Kratos
{
        
  KRATOS_CREATE_VARIABLE(double, AUX_DOUBLE_VAR) //SALVA
  //KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_GRADIENT) //SALVA

KratosSwimmingDEMApplication::KratosSwimmingDEMApplication()
{}

void KratosSwimmingDEMApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosSwimmingDEMApplication... " << std::endl;
                
  KRATOS_REGISTER_VARIABLE(AUX_DOUBLE_VAR) //SALVA
  
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_GRADIENT) //SALVA
          
  //  KRATOS_REGISTER_VARIABLE ( EXPORT_NEIGHBOUR_LIST )


		/* Define In Global variables.cpp
		  */

	Serializer::Register( "VariablesList", mVariablesList );

	
 }

}  // namespace Kratos.


