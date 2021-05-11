//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#include "testing.h"
#include "includes/parallel_environment.h"

namespace Kratos {
namespace Testing {

DataCommunicator& GetDefaultDataCommunicator()
{
    return ParallelEnvironment::GetDefaultDataCommunicator();
}

}
}
