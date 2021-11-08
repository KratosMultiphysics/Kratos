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

void KratosMetisApplication::Register()
{
    KRATOS_INFO("") << "    KRATOS  __  __      _   _\n"
                    << "           |  \\/  | ___| |_(_)___\n"
                    << "           | |\\/| |/ _ \\ __| / __|\n"
                    << "           | |  | |  __/ |_| \\__ \\\n"
                    << "           |_|  |_|\\___|\\__|_|___/\n"
                    << "Initializing KratosMetisApplication..." << std::endl;
}

}  // namespace Kratos.
