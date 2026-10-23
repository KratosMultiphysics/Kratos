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
#include "future/containers/entity_data_container.h"

namespace Kratos::Future
{

EntityDataContainer::EntityDataContainer(std::size_t BufferSize, std::size_t ChunkSize)
    : mDataContainer(ChunkSize),
      mCapacity(ChunkSize),
      mBufferSize(BufferSize)
{
    KRATOS_ERROR_IF(BufferSize == 0) << "Buffer size must be greater than zero." << std::endl;
    KRATOS_ERROR_IF(ChunkSize == 0) << "Chunk size must be greater than zero." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer::SlotType EntityDataContainer::RegisterEntityId(IndexType EntityId)
{
    KRATOS_TRY

    const auto it_slot = mIdToSlot.find(EntityId);
    if (it_slot != mIdToSlot.end()) {
        return it_slot->second; // idempotent
    }

    const SlotType slot = mNextSlot++;
    mIdToSlot.emplace(EntityId, slot);
    EnsureCapacity();
    return slot;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDataContainer::UnregisterEntityId(IndexType EntityId)
{
    KRATOS_TRY

    const auto it_slot = mIdToSlot.find(EntityId);
    KRATOS_ERROR_IF(it_slot == mIdToSlot.end()) << "Entity id " << EntityId << " is not registered in this EntityDataContainer." << std::endl;
    mIdToSlot.erase(it_slot); // the slot becomes a hole; it is not reclaimed (Phase II limitation)

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDataContainer::CloneStepData(StepCategory Category)
{
    mDataContainer.CloneStepData(Category);
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer::SlotType EntityDataContainer::Index(IndexType EntityId) const
{
    KRATOS_TRY

    const auto it_slot = mIdToSlot.find(EntityId);
    KRATOS_ERROR_IF(it_slot == mIdToSlot.end()) << "Entity id " << EntityId << " is not registered in this EntityDataContainer." << std::endl;
    return it_slot->second;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

DataContainer& EntityDataContainer::GetDataContainer()
{
    return mDataContainer;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t EntityDataContainer::GetBufferSize() const
{
    return mBufferSize;
}

/***********************************************************************************/
/***********************************************************************************/

bool EntityDataContainer::HasEntity(IndexType EntityId) const
{
    return mIdToSlot.find(EntityId) != mIdToSlot.end();
}

/***********************************************************************************/
/***********************************************************************************/

bool EntityDataContainer::Has(const VariableData& rVariable) const
{
    return mVariableRecords.find(rVariable.Key()) != mVariableRecords.end();
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t EntityDataContainer::NumberOfEntities() const
{
    return mIdToSlot.size();
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t EntityDataContainer::Capacity() const
{
    return mCapacity;
}

/***********************************************************************************/
/***********************************************************************************/

std::string EntityDataContainer::Info() const
{
    return "EntityDataContainer";
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDataContainer::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDataContainer::PrintData(std::ostream& rOStream) const
{
    rOStream << "Registered entities: " << NumberOfEntities() << std::endl;
    rOStream << "Slot capacity      : " << Capacity() << std::endl;
    rOStream << "Buffer size        : " << GetBufferSize() << std::endl;
    rOStream << "Added variables    : " << mVariableRecords.size() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDataContainer::EnsureCapacity()
{
    KRATOS_TRY

    if (mNextSlot <= mCapacity) {
        return;
    }

    mCapacity = std::max(2 * mCapacity, static_cast<std::size_t>(mNextSlot));
    mDataContainer.Resize(mCapacity);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const EntityDataContainer::VariableRecord& EntityDataContainer::GetVariableRecord(const VariableData& rVariable, int StepBeforeCurrent) const
{
    KRATOS_TRY

    const auto it_record = mVariableRecords.find(rVariable.Key());
    KRATOS_ERROR_IF(it_record == mVariableRecords.end()) << "Variable " << rVariable.Name() << " has not been added to this EntityDataContainer." << std::endl;
    KRATOS_ERROR_IF(StepBeforeCurrent != 0 && it_record->second.Category == StepCategory::AnyStep)
        << "Variable " << rVariable.Name() << " is not historical: requesting " << StepBeforeCurrent << " steps before the current one is not possible." << std::endl;
    return it_record->second;

    KRATOS_CATCH("")
}

} // namespace Kratos::Future
