//  ______  _                _            ______        ______                           _ _                 _
// (_____ \| |              (_)          |  ___ \      |  ___ \         /\              | (_)           _   (_)
//  _____) ) | _  _   _  ___ _  ____  ___| |   | | ____| | _ | | ___   /  \  ____  ____ | |_  ____ ____| |_  _  ___  ____
// |  ____/| || \| | | |/___) |/ ___)/___) |   | |/ _  ) || || |/ _ \ / /\ \|  _ \|  _ \| | |/ ___) _  |  _)| |/ _ \|  _ \
// | |     | | | | |_| |___ | ( (___|___ | |   | ( (/ /| || || | |_| | |__| | | | | | | | | ( (__( ( | | |__| | |_| | | | |
// |_|     |_| |_|\__  (___/|_|\____|___/|_|   |_|\____)_||_||_|\___/|______| ||_/| ||_/|_|_|\____)_||_|\___)_|\___/|_| |_|
//               (____/                                                     |_|   |_|
//
//  License:         BSD License
//                   license: PhysicsNeMoApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <sstream>

// External includes

// Project includes
#include "physics_nemo_application.h"
#include "physics_nemo_application_variables.h"

namespace Kratos {

KratosPhysicsNeMoApplication::KratosPhysicsNeMoApplication()
    : KratosApplication("PhysicsNeMoApplication")
{
}

void KratosPhysicsNeMoApplication::Register()
{
    std::stringstream banner;

    banner << R"( ______  _                _            ______        ______                           _ _                 _            )" << "\n"
           << R"((_____ \| |              (_)          |  ___ \      |  ___ \         /\              | (_)           _   (_)            )" << "\n"
           << R"( _____) ) | _  _   _  ___ _  ____  ___| |   | | ____| | _ | | ___   /  \  ____  ____ | |_  ____ ____| |_  _  ___  ____  )" << "\n"
           << R"(|  ____/| || \| | | |/___) |/ ___)/___) |   | |/ _  ) || || |/ _ \ / /\ \|  _ \|  _ \| | |/ ___) _  |  _)| |/ _ \|  _ \)" << "\n"
           << R"(| |     | | | | |_| |___ | ( (___|___ | |   | ( (/ /| || || | |_| | |__| | | | | | | | | ( (__( ( | | |__| | |_| | | | |)" << "\n"
           << R"(|_|     |_| |_|\__  (___/|_|\____|___/|_|   |_|\____)_||_||_|\___/|______| ||_/| ||_/|_|_|\____)_||_|\___)_|\___/|_| |_|)" << "\n"
           << R"(              (____/                                                     |_|   |_|                )" << "\n"
           << "Initializing KratosPhysicsNeMoApplication..." << std::endl;

    KRATOS_INFO("") << banner.str();
}

} // namespace Kratos
