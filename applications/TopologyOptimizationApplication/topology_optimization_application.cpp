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
#include "topology_optimization_application.h"
#include "topology_optimization_application_variables.h"


namespace Kratos {

KratosTopologyOptimizationApplication::KratosTopologyOptimizationApplication():
    KratosApplication("TopologyOptimizationApplication")
    {}

void KratosTopologyOptimizationApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosTopologyOptimizationApplication..." << std::endl;
}
}  // namespace Kratos.
