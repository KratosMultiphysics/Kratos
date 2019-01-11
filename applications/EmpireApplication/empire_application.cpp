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
#include "empire_application.h"
#include "includes/variables.h"


namespace Kratos
{

KratosEmpireApplication::KratosEmpireApplication() : KratosApplication("EmpireApplication") {}

void KratosEmpireApplication::Register()
{
   // calling base class register to register Kratos components
   KratosApplication::Register();

   KRATOS_INFO("") << "    KRATOS  _____ __  __ ____ ___ ____  _____\n"
                   << "           | ____|  \\/  |  _ \\_ _|  _ \\| ____|\n"
                   << "           |  _| | |\\/| | |_) | || |_) |  _|\n"
                   << "           | |___| |  | |  __/| ||  _ <| |___\n"
                   << "           |_____|_|  |_|_|  |___|_| \\_\\_____|\n"
                   << "Initializing KratosEmpireApplication..." << std::endl;

}

}  // namespace Kratos.

