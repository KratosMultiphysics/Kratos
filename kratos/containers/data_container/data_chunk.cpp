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
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "containers/data_container/data_chunk.h"

namespace Kratos
{

void DataChunkBase::ResizeData(std::size_t NumberOfEntities)
{
    this->ResizeDataContainer(NumberOfEntities);
    mNumberOfEntitiesPerStep = NumberOfEntities;
}

} // namespace Kratos
