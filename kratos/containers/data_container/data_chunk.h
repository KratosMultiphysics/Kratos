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

#pragma once

// System includes
#include <algorithm>

// External includes

// Project includes
#include "containers/data_container/data_value_policy.h"
#include "containers/data_container/data_history_policy.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DataChunkBase
 * @ingroup KratosCore
 * @brief Type-erased base class of the storage chunks held by a DataContainer.
 * @details Provides the type-independent part of a chunk: the number of entities stored per step and the virtual interface used by DataContainer to work with chunks without knowing their value type (query the associated variable and policies, create sibling chunks, clone step data and resize the storage).
 *
 * Chunks are usually constructed and owned (as std::shared_ptr) by DataContainer::Add, not directly. A chunk with zero entities is valid: chunks governed by a sparse value policy start empty and grow via DataContainer::UpdateSparseStorage / DataContainer::AddToSparseStorage.
 * @note This class provides no thread-safety guarantees; concurrent access must be synchronized externally (DataContainer only locks chunk creation in Add).
 * @see DataChunk, DataContainer, DataValuePolicyBase, DataHistoryPolicyBase
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) DataChunkBase
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param NumberOfEntities Number of entities per step in the chunk. May be zero for (initially empty) sparse chunks.
     */
    DataChunkBase(std::size_t NumberOfEntities)
        : mNumberOfEntitiesPerStep(NumberOfEntities)
    {
    }

    /// Destructor.
    virtual ~DataChunkBase() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a new chunk of the same dynamic type, variable and policies with the given number of entities.
     * @param NumberOfEntities Number of entities per step of the new chunk.
     * @return A shared pointer to the newly created chunk (its data is zero-initialized and its history ring starts at slot 0, see HistoricalDataPolicy::Clone).
     */
    virtual std::shared_ptr<DataChunkBase> CreateNew(std::size_t NumberOfEntities) const = 0;

    /**
     * @brief Allocate the raw storage for the given number of entities.
     * @details Any previously held data is released (without copying); the new storage is uninitialized. Allocates NumberOfEntities times the total number of steps values.
     * @param NumberOfEntities Number of entities per step to allocate for. May be zero (no storage is allocated).
     */
    virtual void AllocateData(std::size_t NumberOfEntities) = 0;

    /**
     * @brief Clone the data of the current step onto the next step slot of the given category.
     * @details Delegates the ring-buffer advance to the history policy: data is copied only when the policy's category matches and the step index actually changes; otherwise this is a no-op.
     * @param Category The step category being cloned.
     */
    virtual void CloneStepData(StepCategory Category) = 0;

    /**
     * @brief Reallocate the storage preserving existing data.
     * @details For each stored step the first min(old, new) entity values are preserved; newly added entries are zero-initialized with the value policy zero. Resizing to zero releases the storage. This reallocates and copies (O(entities * steps)).
     * @param NumberOfEntities The new number of entities per step.
     */
    void ResizeData(std::size_t NumberOfEntities);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the variable associated with the data chunk.
     * @return Reference to the variable data.
     */
    virtual const VariableData& GetVariableData() const = 0;

    /**
     * @brief Get the value policy associated with the data chunk.
     * @return Reference to the value policy.
     */
    virtual const DataValuePolicyBase& GetValuePolicy() const = 0;

    /**
     * @brief Get the history policy associated with the data chunk.
     * @return Reference to the history policy.
     */
    virtual const DataHistoryPolicyBase& GetHistoryPolicy() const = 0;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the number of entities stored per step.
     * @return The number of entities per step.
     */
    std::size_t NumberOfEntitiesPerStep() const
    {
        return mNumberOfEntitiesPerStep;
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief Reallocate the storage for the new number of entities preserving existing data.
     * @details Implemented by DataChunk; called by ResizeData, which afterwards updates the stored number of entities.
     * @param NumberOfEntities The new number of entities per step.
     */
    virtual void ResizeDataContainer(std::size_t NumberOfEntities) = 0;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::size_t mNumberOfEntitiesPerStep; /// Number of entities per step in the chunk. This is used to determine the size of the data

    ///@}

}; // Class DataChunkBase

/**
 * @class DataChunk
 * @ingroup KratosCore
 * @brief Contiguous storage for the values of one variable over the entities of a DataContainer.
 * @details The chunk owns a single raw array of NumberOfEntitiesPerStep() times GetTotalNumberOfSteps() values of type TValueType, laid out step-by-step (all entities of step slot 0, then all entities of step slot 1, ...). Ownership is manual (new[] / delete[]) to keep the storage a plain contiguous block, mirroring how other performance-critical Kratos containers manage their memory; the chunk is therefore non-copyable (it is shared via std::shared_ptr from the DataContainer instead).
 *
 * On construction all values are initialized with the value policy zero. The value policy and history policy passed in are cloned and owned by the chunk; the variable is held by reference and must outlive the chunk (DataContainer passes the Registry-owned canonical variable).
 * @note This class provides no thread-safety guarantees.
 * @see DataChunkBase, DataContainer, DataValuePolicy, DataHistoryPolicyBase
 * @tparam TValueType Type of the stored values.
 * @author Pooyan Dadvand
 */
template <typename TValueType>
class DataChunk : public DataChunkBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Type of the stored values
    using ValueType = TValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rVariable The variable associated with the data chunk. Held by reference; it must outlive the chunk.
     * @param rValuePolicy The value policy; it is cloned into the chunk.
     * @param rHistoryPolicy The history policy; it is cloned into the chunk.
     * @param NumberOfEntities Number of entities per step. May be zero for sparse chunks.
     */
    DataChunk(
        const VariableData& rVariable,
        const DataValuePolicy<TValueType>& rValuePolicy,
        const DataHistoryPolicyBase& rHistoryPolicy,
        std::size_t NumberOfEntities)
        : DataChunkBase(NumberOfEntities),
          mData(nullptr),
          mpValuePolicy(rValuePolicy.Clone()),
          mpHistoryPolicy(rHistoryPolicy.Clone()),
          mrVariable(rVariable)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpHistoryPolicy && mpHistoryPolicy->GetTotalNumberOfSteps() > 0) << "Number of buffers must be greater than zero." << std::endl;

        AllocateData(NumberOfEntities);
        if (mData != nullptr) {
            std::fill(mData, mData + NumberOfEntities * mpHistoryPolicy->GetTotalNumberOfSteps(), rValuePolicy.Zero());
        }
    }

    /**
     * @brief Constructor taking ownership of already cloned policies.
     * @param rVariable The variable associated with the data chunk. Held by reference; it must outlive the chunk.
     * @param pValuePolicy The value policy; ownership is transferred to the chunk.
     * @param pHistoryPolicy The history policy; ownership is transferred to the chunk.
     * @param NumberOfEntities Number of entities per step. May be zero for sparse chunks.
     */
    template <typename TValuePolicyType>
    DataChunk(
        const VariableData& rVariable,
        std::unique_ptr<TValuePolicyType> pValuePolicy,
        std::unique_ptr<DataHistoryPolicyBase> pHistoryPolicy,
        std::size_t NumberOfEntities)
        : DataChunkBase(NumberOfEntities),
          mData(nullptr),
          mpValuePolicy(std::move(pValuePolicy)),
          mpHistoryPolicy(std::move(pHistoryPolicy)),
          mrVariable(rVariable)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpHistoryPolicy && mpHistoryPolicy->GetTotalNumberOfSteps() > 0) << "Number of buffers must be greater than zero." << std::endl;

        AllocateData(NumberOfEntities);
        if (mData != nullptr) {
            std::fill(mData, mData + NumberOfEntities * mpHistoryPolicy->GetTotalNumberOfSteps(), mpValuePolicy->Zero());
        }
    }

    /// Copy constructor (deleted: the chunk owns raw storage and is shared via std::shared_ptr instead).
    DataChunk(const DataChunk&) = delete;

    /// Destructor.
    ~DataChunk() override
    {
        delete[] mData;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator (deleted: the chunk owns raw storage and is shared via std::shared_ptr instead).
    DataChunk& operator=(const DataChunk&) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a new DataChunk of the same type, variable and policies to store NumberOfEntities objects.
     * @param NumberOfEntities Number of entities per step of the new chunk.
     * @return A shared pointer to the newly created, zero-initialized chunk. Its history ring starts at slot 0 (the history policy is cloned, see HistoricalDataPolicy::Clone).
     */
    std::shared_ptr<DataChunkBase> CreateNew(std::size_t NumberOfEntities) const override
    {
        return std::make_shared<DataChunk<TValueType>>(mrVariable, mpValuePolicy->Clone(), mpHistoryPolicy->Clone(), NumberOfEntities);
    }

    /**
     * @brief Allocate the raw storage for the given number of entities.
     * @details Releases any previously held data (without copying); the new storage is uninitialized (values are default-constructed by new[]). Allocates NumberOfEntities times the total number of steps values; with zero entities no storage is allocated.
     * @param NumberOfEntities Number of entities per step to allocate for.
     */
    void AllocateData(std::size_t NumberOfEntities) override
    {
        delete[] mData;
        mData = nullptr;
        const std::size_t total_size = NumberOfEntities * mpHistoryPolicy->GetTotalNumberOfSteps();
        if (total_size > 0) {
            mData = new ValueType[total_size];
        }
    }

    /**
     * @brief Clone the data of the current step onto the next step slot of the given category.
     * @details The history policy advances its ring only when the category matches its own; when the step index changes, the values of the previous current step are copied onto the new current step (so the new step starts as a copy of the old one). Otherwise this is a no-op.
     * @param Category The step category being cloned.
     */
    void CloneStepData(StepCategory Category) override
    {
        const std::size_t current_step_index = mpHistoryPolicy->GetLatestStepIndex();
        const std::size_t next_step_index = mpHistoryPolicy->CloneStepIndex(Category);

        if (next_step_index != current_step_index) {
            const std::size_t offset_current = current_step_index * NumberOfEntitiesPerStep();
            const std::size_t offset_next = next_step_index * NumberOfEntitiesPerStep();
            std::copy(mData + offset_current, mData + offset_current + NumberOfEntitiesPerStep(),
                      mData + offset_next);
        }
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the variable associated with the data chunk.
     * @return Reference to the variable data.
     */
    const VariableData& GetVariableData() const override
    {
        return mrVariable;
    }

    /**
     * @brief Get the value policy associated with the data chunk.
     * @return Reference to the value policy.
     */
    const DataValuePolicy<TValueType>& GetValuePolicy() const override
    {
        return *mpValuePolicy;
    }

    /**
     * @brief Get the history policy associated with the data chunk.
     * @return Reference to the history policy.
     */
    const DataHistoryPolicyBase& GetHistoryPolicy() const override
    {
        return *mpHistoryPolicy;
    }

    /**
     * @brief Get the raw data pointer (all steps).
     * @return Pointer to the first value of step slot 0 (nullptr for empty chunks).
     */
    ValueType* GetData()
    {
        return mData;
    }

    /**
     * @brief Get the raw data pointer (all steps, const version).
     * @return Pointer to the first value of step slot 0 (nullptr for empty chunks).
     */
    const ValueType* GetData() const
    {
        return mData;
    }

    /**
     * @brief Get the iterator to the data of the given step slot.
     * @param StepIndex The step slot index (see DataHistoryPolicyBase::GetStepIndex).
     * @return Pointer to the first value of the given step slot.
     */
    ValueType* GetDataIterator(std::size_t StepIndex)
    {
        KRATOS_DEBUG_ERROR_IF(StepIndex >= mpHistoryPolicy->GetTotalNumberOfSteps()) << "Step index out of bounds." << std::endl;
        const std::size_t offset = StepIndex * NumberOfEntitiesPerStep();
        return mData + offset;
    }

    /**
     * @brief Get the pointer to the data of the given step slot.
     * @details Same as GetDataIterator (kept for interface compatibility).
     * @param StepIndex The step slot index (see DataHistoryPolicyBase::GetStepIndex).
     * @return Pointer to the first value of the given step slot.
     */
    ValueType* GetDataPointer(std::size_t StepIndex)
    {
        KRATOS_DEBUG_ERROR_IF(StepIndex >= mpHistoryPolicy->GetTotalNumberOfSteps()) << "Step index out of bounds." << std::endl;
        const std::size_t offset = StepIndex * NumberOfEntitiesPerStep();
        return mData + offset;
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief Reallocate the storage preserving existing data.
     * @details For each stored step slot the first min(old, new) entity values are copied into the new storage and newly added entries are initialized with the value policy zero, so per-step data is preserved independently of the history depth. Resizing to zero releases the storage. Reallocates and copies (O(entities * steps)).
     * @param NewNumberOfEntities The new number of entities per step.
     */
    void ResizeDataContainer(std::size_t NewNumberOfEntities) override
    {
        if (NewNumberOfEntities == 0) {
            delete[] mData;
            mData = nullptr;
            return;
        }

        const std::size_t num_steps = mpHistoryPolicy->GetTotalNumberOfSteps();
        const std::size_t current_num_entities = NumberOfEntitiesPerStep();
        // Per-step number of preserved entries (the old data may be empty)
        const std::size_t copy_entities = (mData != nullptr) ? std::min(current_num_entities, NewNumberOfEntities) : 0;

        ValueType* p_data = mData;
        mData = new ValueType[NewNumberOfEntities * num_steps];

        for (std::size_t step = 0; step < num_steps; step++) {
            if (copy_entities > 0) {
                std::copy(
                    p_data + step * current_num_entities,
                    p_data + step * current_num_entities + copy_entities,
                    mData + step * NewNumberOfEntities);
            }

            std::fill(
                mData + step * NewNumberOfEntities + copy_entities,
                mData + (step + 1) * NewNumberOfEntities,
                mpValuePolicy->Zero());
        }

        delete[] p_data;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    TValueType* mData = nullptr;                              /// The data stored in the chunk as a raw owning pointer

    std::unique_ptr<DataValuePolicy<TValueType>> mpValuePolicy; /// The value policy associated with the data chunk. This is used to determine the type of the data to store

    std::unique_ptr<DataHistoryPolicyBase> mpHistoryPolicy;   /// Storing the history policy as a unique pointer allows for polymorphic behavior and memory management

    const VariableData& mrVariable;                           /// The variable associated with the data chunk

    ///@}

}; // Class DataChunk

///@}

///@} addtogroup block

} // namespace Kratos
