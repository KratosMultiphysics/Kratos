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
#include "metis_application.h"

namespace Kratos {

KratosMetisApplication::KratosMetisApplication()
    : KratosApplication("MetisApplication") {}

void KratosMetisApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing Kratos MetisApplication... " << std::endl;
}

}  // namespace Kratos.
