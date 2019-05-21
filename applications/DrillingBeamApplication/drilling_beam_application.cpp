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
#include "drilling_beam_application.h"

namespace Kratos
{

	//KRATOS_CREATE_VARIABLE(double, NODAL_AREA)


 	KratosDrillingBeamApplication::KratosDrillingBeamApplication():
	 KratosApplication("DrillingBeamApplication")
 	{}

 	void KratosDrillingBeamApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
        KRATOS_INFO("") << "Initializing KratosDrillingBeamApplication... " << std::endl;

// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 	}

}  // namespace Kratos.


