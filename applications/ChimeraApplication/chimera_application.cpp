//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//

// System includes

// External includes

// Project includes
#include "chimera_application.h"
#include "chimera_application_variables.h"

namespace Kratos
{

KratosChimeraApplication::KratosChimeraApplication():KratosApplication("ChimeraApplication")
{}

void KratosChimeraApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") <<"Kratos       ________  ________  _____________  ___   \n"
                    <<"            / ____/ / / /  _/  |/  / ____/ __ \\/   | \n"
                    <<"           / /   / /_/ // // /|_/ / __/ / /_/ / /| |  \n"
                    <<"          / /___/ __  // // /  / / /___/ _, _/ ___ |  \n"
                    <<"          \\____/_/ /_/___/_/  /_/_____/_/ |_/_/  |_| \n"
                    <<"                                                     Application initializing ... "<<std::endl;
    // Flag for distinguishing b/w velocity and pressure constraints.
    KRATOS_REGISTER_FLAG(FS_CHIMERA_VEL_CONSTRAINT);
    KRATOS_REGISTER_FLAG(FS_CHIMERA_PRE_CONSTRAINT);
    KRATOS_REGISTER_FLAG(CHIMERA_INTERNAL_BOUNDARY);
}

} // namespace Kratos.
