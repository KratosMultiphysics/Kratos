//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"

namespace Kratos {

DataCommunicator& DataCommunicator::GetDefault()
{
    return ParallelEnvironment::GetDefaultDataCommunicator();
}

}