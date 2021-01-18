//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    janosch
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "externalsolvers_application.h"


namespace Kratos
{

void KratosExternalSolversApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosExternalSolversApplication..." << std::endl;

    ExternalSolversApplicationRegisterLinearSolvers();

}


}  // namespace Kratos.


