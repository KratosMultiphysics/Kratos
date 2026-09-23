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
#include "containers/data_container/data_container.h"

namespace Kratos
{

DataContainer::DataContainer(std::size_t ChunkSize)
    : mChunkSize(ChunkSize)
{
    KRATOS_DEBUG_ERROR_IF(mChunkSize == 0) << "Chunk size must be greater than zero." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

DataContainer::DataContainer(const DataContainer& rOther)
    : mData(rOther.mData),
      mChunkSize(rOther.mChunkSize)
{
}

/***********************************************************************************/
/***********************************************************************************/

DataContainer::DataContainer(DataContainer&& rOther) noexcept
    : mData(std::move(rOther.mData)),
      mChunkSize(rOther.mChunkSize)
{
}

/***********************************************************************************/
/***********************************************************************************/

DataContainer& DataContainer::operator=(const DataContainer& rOther)
{
    mData = rOther.mData;
    mChunkSize = rOther.mChunkSize;
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

DataContainer& DataContainer::operator=(DataContainer&& rOther) noexcept
{
    mData = std::move(rOther.mData);
    mChunkSize = rOther.mChunkSize;
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::Initialize(const DataContainer& rOther, std::size_t ChunkSize)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(ChunkSize == 0) << "Chunk size must be greater than zero." << std::endl;
    mChunkSize = ChunkSize;

    // Copy chunks from the other data container
    for (const auto& p_chunk : rOther.mData) {
        if (p_chunk->GetValuePolicy().IsSparse()) {
            mData.push_back(p_chunk->CreateNew(0));
        } else {
            mData.push_back(p_chunk->CreateNew(ChunkSize));
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::UpdateSparseStorage(const DataAccessor<int>& rIndexAccessor)
{
    KRATOS_TRY

    for (auto& p_chunk : mData) {
        if (p_chunk->GetValuePolicy().IsSparse()) {
            auto index_accessor = p_chunk->GetValuePolicy().GetSparseIndexAccessor();
            if (index_accessor == rIndexAccessor) {
                auto index_span = GetDataSpan(index_accessor);
                std::size_t sparse_count = 0;
                for (const auto index : index_span) {
                    if (index > -1) {
                        ++sparse_count;
                    }
                }

                auto p_new_chunk = p_chunk->CreateNew(sparse_count);
                p_chunk.swap(p_new_chunk);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::AddToSparseStorage(const DataAccessor<int>& rIndexAccessor, const std::vector<int>& rEntityIndices)
{
    KRATOS_TRY

    if (rEntityIndices.size() == 0) return;

    int new_size = -1;
    for (auto& p_chunk : mData) {
        if (p_chunk->GetValuePolicy().IsSparse()) {
            auto index_accessor = p_chunk->GetValuePolicy().GetSparseIndexAccessor();
            if (index_accessor == rIndexAccessor) {
                if (new_size < 0) {
                    new_size = static_cast<int>(p_chunk->NumberOfEntitiesPerStep());
                    auto index_span = GetDataSpan(index_accessor);
                    for (auto entity_id : rEntityIndices) {
                        if (index_span[entity_id] == -1) {
                            index_span[entity_id] = new_size++;
                        }
                    }
                }

                p_chunk->ResizeData(new_size);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::Resize(std::size_t NewNumberOfEntities)
{
    KRATOS_TRY

    for (auto& p_chunk : mData) {
        if (p_chunk->GetValuePolicy().IsSparse()) continue; // sparse sizing is governed by UpdateSparseStorage/AddToSparseStorage
        if (p_chunk->NumberOfEntitiesPerStep() == NewNumberOfEntities) continue;
        p_chunk->ResizeData(NewNumberOfEntities);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::CloneStepData(StepCategory rStepCategory)
{
    for (auto& p_chunk : mData) {
        p_chunk->CloneStepData(rStepCategory);
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::string DataContainer::Info() const
{
    return "DataContainer";
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void DataContainer::PrintData(std::ostream& rOStream) const
{
    rOStream << "Number of chunks: " << mData.size() << std::endl;
    rOStream << "Chunk size      : " << mChunkSize << std::endl;
    for (const auto& p_chunk : mData) {
        rOStream << "    " << p_chunk->GetVariableData().Name()
                 << " (" << p_chunk->NumberOfEntitiesPerStep() << " entities per step, "
                 << p_chunk->GetHistoryPolicy().GetTotalNumberOfSteps() << " steps)" << std::endl;
    }
}

} // namespace Kratos
