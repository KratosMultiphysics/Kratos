//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
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
     KRATOS_INFO("") << "Initializing KratosNeuralNetworkApplication..." << std::endl;


}
}  // namespace Kratos.
