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

void KratosFSIApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosFSIApplication... " << std::endl;

}

}  // namespace Kratos.
