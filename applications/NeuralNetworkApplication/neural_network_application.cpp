//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rishith Ellath Meethal (https://github.com/rishithellathmeethal)
//                   Daniel Andrés Arcones https://github.com/danielandresarcones
//

// System includes

// External includes

// Project includes
#include "neural_network_application.h"
#include "neural_network_application_variables.h"

namespace Kratos {

KratosNeuralNetworkApplication::KratosNeuralNetworkApplication():
    KratosApplication("NeuralNetworkApplication")
    {}

void KratosNeuralNetworkApplication::Register()
{
     KRATOS_INFO("") 
     << " KRATOS     _   __                      ___   __     __                      __   ___                ___            __  _            \n"
     << "           / | / /__  __  ___________ _/ / | / /__  / /__      ______  _____/ /__/   |  ____  ____  / (_)________ _/ /_(_)___  ____  \n"
     << "          /  |/ / _ \/ / / / ___/ __ `/ /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/ /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \ \n"
     << "         / /|  /  __/ /_/ / /  / /_/ / / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,< / ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / /_/ / / / / \n"
     << "        /_/ |_/\___/\__,_/_/   \__,_/_/_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_/_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/  \n"
     << "                                                                                    /_/   /_/                                        \n"
     << "Initializing KratosNeuralNetworkApplication..." << std::endl;

}
}  // namespace Kratos.
