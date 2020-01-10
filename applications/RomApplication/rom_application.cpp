//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:      BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors: Riccardo Rossi
//                Raul Bravo
//


// System includes


// External includes


// Project includes
#include "rom_application.h"
#include "rom_application_variables.h"


namespace Kratos {

KratosRomApplication::KratosRomApplication():
    KratosApplication("RomApplication")
    {}

void KratosRomApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosRomApplication..." << std::endl;
    KRATOS_REGISTER_VARIABLE( AUX_ID )
    KRATOS_REGISTER_VARIABLE( ROM_BASIS )
}
}  // namespace Kratos.
