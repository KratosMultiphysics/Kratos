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

#ifndef KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED

#include <string>

#include "includes/data_communicator.h"

namespace Kratos
{

namespace DataCommunicatorFactory
{

const DataCommunicator& DuplicateAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    const std::string& rNewCommunicatorName);

const DataCommunicator& SplitAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    int Color,
    int Key,
    const std::string& rNewCommunicatorName);

const DataCommunicator& CreateFromRanksAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    const std::vector<int>& rRanks,
    const std::string& rNewCommunicatorName);

const DataCommunicator& CreateUnionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName);

const DataCommunicator& CreateIntersectionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName);

}

}

#endif // KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED