//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "kahip_application.h"

namespace Kratos {

KratosKaHIPApplication::KratosKaHIPApplication()
    : KratosApplication("KaHIPApplication") {}

void KratosKaHIPApplication::Register()
{
    KRATOS_INFO("") << "    KRATOS   _  __    _  _ ___ ____\n"
                    << "            | |/ /__ _| || |_ _|  _ \\\n"
                    << "            | ' // _` | __ || || |_) |\n"
                    << "            | . \\ (_| | || || ||  __/\n"
                    << "            |_|\\_\\__,_|_||_|___|_|\n"
                    << "Initializing KratosKaHIPApplication..." << std::endl;
}

}  // namespace Kratos
