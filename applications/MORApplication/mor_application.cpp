//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "mor_application.h"
#include "mor_application_variables.h"


namespace Kratos {

KratosMORApplication::KratosMORApplication():
    KratosApplication("MORApplication")
    {}

void KratosMORApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosMORApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( FREQUENCY )

}
}  // namespace Kratos.
