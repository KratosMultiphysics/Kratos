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

// External includes

// Project includes
#include "includes/lock_object.h"
#include "includes/registry.h"
#include "containers/variable.h"
#include "containers/data_container/data_accessor.h"
#include "containers/data_container/data_chunk.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DataContainer
 * @ingroup KratosCore
 * @brief Chunked, type-erased, variable-keyed storage of entity data.
 * @details The container stores the data associated with each added Variable in one DataChunk: a contiguous block holding the values of that variable for all entities (struct-of-arrays layout) and, depending on the chunk's history policy, for several step slots. How values are stored is governed by a DataValuePolicy (dense, layered or sparse) and how step history is kept by a DataHistoryPolicy, both passed to Add.
 *
 * Add returns a DataAccessor caching the chunk index, giving O(1) access to the data span later on (GetDataSpan). Variables passed to Add are canonicalized through the Kratos Registry ("variables.all.<name>"): the chunks and accessors reference the Registry-owned variable, so the stored data remains valid and retrievable even if the Variable object passed by the caller goes out of scope. Consequently, only registered variables can be added (see KRATOS_REGISTER_VARIABLE / Variable::Register); Add does not register variables itself, it only canonicalizes already-registered ones.
 *
 * Dense chunks are created with mChunkSize entities; chunks with a sparse value policy start with zero entities and grow via UpdateSparseStorage / AddToSparseStorage.
 *
 * Copying a DataContainer is shallow: the copy shares the chunks (std::shared_ptr) with the original. Thread-safety: only Add is internally locked (chunk creation); all other operations require external synchronization.
 * @see DataChunk, DataAccessor, DataValuePolicy, SparseDataValuePolicy, DataHistoryPolicyBase
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) DataContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataContainer
    KRATOS_CLASS_POINTER_DEFINITION(DataContainer);

    /// Accessor type alias for the stored data of a given value type
    template<typename TValueType>
    using DataAccessor = Kratos::DataAccessor<TValueType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param ChunkSize Number of entities each dense chunk stores per step. Must be greater than zero.
     */
    DataContainer(std::size_t ChunkSize = 256);

    /// Destructor.
    ~DataContainer() = default;

    /**
     * @brief Copy constructor.
     * @details The copy is shallow: it shares the chunks with rOther. The lock is default-constructed (locks are not copyable).
     * @param rOther The container to copy from.
     */
    DataContainer(const DataContainer& rOther);

    /**
     * @brief Move constructor.
     * @details The lock is default-constructed (locks are not movable); rOther is left without chunks.
     * @param rOther The container to move from.
     */
    DataContainer(DataContainer&& rOther) noexcept;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Copy assignment operator.
     * @details The copy is shallow: this container shares the chunks with rOther afterwards. The lock of each container is untouched.
     * @param rOther The container to copy from.
     * @return Reference to this container.
     */
    DataContainer& operator=(const DataContainer& rOther);

    /**
     * @brief Move assignment operator.
     * @param rOther The container to move from; it is left without chunks.
     * @return Reference to this container.
     */
    DataContainer& operator=(DataContainer&& rOther) noexcept;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize this container with fresh chunks mirroring another container.
     * @details For every chunk of rOther a new sibling chunk (same variable and policies) is appended to this container: dense chunks are created with the given chunk size, sparse chunks start empty. The data itself is not copied (all new chunks are zero-initialized) and the given chunk size is stored for future Add calls.
     * @param rOther The container whose chunk structure is mirrored.
     * @param ChunkSize Number of entities per step of the new dense chunks. Must be greater than zero.
     */
    void Initialize(const DataContainer& rOther, std::size_t ChunkSize);

    /**
     * @brief Add a data chunk for the given variable and policies.
     * @details The variable is canonicalized through the Kratos Registry ("variables.all.<name>"), so it must be registered beforehand (see KRATOS_REGISTER_VARIABLE / Variable::Register); the chunk and the returned accessor reference the Registry-owned variable, making the stored data independent of the lifetime of rVariable. If the variable was already added with the same value policy (IsSame), the accessor of the existing chunk is returned; adding it with a different policy is an error. New dense chunks are created with the container chunk size, sparse chunks start with zero entities.
     * @note This method is internally locked and safe to call from multiple threads.
     * @tparam TVariableValueType The value type of the variable.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable for which the data chunk is to be added.
     * @param rValuePolicy The value policy governing the stored values (cloned into the chunk).
     * @param rHistoryPolicy The history policy governing the step slots (cloned into the chunk).
     * @return An accessor to the data of the added (or already existing) chunk.
     */
    template<typename TVariableValueType, typename TValuePolicyType>
    DataAccessor<typename TValuePolicyType::ValueType> Add(
        const Variable<TVariableValueType>& rVariable,
        const TValuePolicyType& rValuePolicy,
        const DataHistoryPolicyBase& rHistoryPolicy = NonHistoricalDataPolicy())
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(rValuePolicy.IsCompatible(rVariable)) << rVariable.Name() << " is not compatible with the given value policy." << std::endl;

        // Guard for mutual exclusion if this function is called from multiple threads
        std::lock_guard<LockObject> lock(mMutex);

        // Always use the reference variable from the registry so that we do not need to care if rVariable goes out of scope
        using VariableType = Variable<TVariableValueType>;

        const std::string variable_path = "variables.all." + rVariable.Name();
        KRATOS_ERROR_IF_NOT(Registry::HasItem(variable_path)) << "Variable " << rVariable.Name()
            << " is not registered in the Registry. DataContainer only accepts registered variables"
            << " (see KRATOS_REGISTER_VARIABLE or Variable::Register)." << std::endl;
        const VariableType& r_reference_variable = Registry::GetValue<VariableType>(variable_path);

        // Check if the variable already exists in the container
        for (auto i_chunk = mData.begin(); i_chunk != mData.end(); ++i_chunk) {
            if ((*i_chunk)->GetVariableData() == r_reference_variable) {
                // If the variable exists, check if the policy matches
                if ((*i_chunk)->GetValuePolicy().IsSame(rValuePolicy)) {
                    return DataAccessor<typename TValuePolicyType::ValueType>(r_reference_variable, std::distance(mData.begin(), i_chunk),
                                                                              rHistoryPolicy.GetStepCategory());
                } else {
                    KRATOS_ERROR << "Variable " << r_reference_variable.Name() << " exists with a different policy." << std::endl;
                }
            }
        }

        // Variable does not exist, create a new data chunk and add it to the container
        const std::size_t initial_chunk_size = rValuePolicy.IsSparse() ? 0 : mChunkSize;
        auto p_new_chunk = std::make_shared<DataChunk<typename TValuePolicyType::ValueType>>(r_reference_variable, rValuePolicy, rHistoryPolicy, initial_chunk_size);
        mData.push_back(p_new_chunk);
        return DataAccessor<typename TValuePolicyType::ValueType>(r_reference_variable, mData.size() - 1,
                                                                  rHistoryPolicy.GetStepCategory());

        KRATOS_CATCH("")
    }

    /**
     * @brief Get an accessor for existing data with matching variable and policy types.
     * @details Policy details do not need to match, only the types (IsSameType): e.g. an accessor for data added with HistoricalDataPolicy(TimeStep, 3) can be obtained with HistoricalDataPolicy(TimeStep, 2). Raises an error if no matching chunk exists.
     * @tparam TVariableType The variable type.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable of the requested data.
     * @param rValuePolicy A value policy of the type the data was added with.
     * @param rHistoryPolicy A history policy of the type the data was added with.
     * @return An accessor to the data of the matching chunk.
     */
    template<typename TVariableType, typename TValuePolicyType>
    DataAccessor<typename TValuePolicyType::ValueType> GetAccessor(
        const TVariableType& rVariable,
        const TValuePolicyType& rValuePolicy,
        const DataHistoryPolicyBase& rHistoryPolicy = NonHistoricalDataPolicy())
    {
        for (auto i_chunk = mData.begin(); i_chunk != mData.end(); ++i_chunk) {
            if ((*i_chunk)->GetVariableData() == rVariable &&
                (*i_chunk)->GetValuePolicy().IsSameType(rValuePolicy)) {
                auto const& r_stored_history_policy = (*i_chunk)->GetHistoryPolicy();
                if (r_stored_history_policy.IsSameType(rHistoryPolicy)) {
                    return DataAccessor<typename TValuePolicyType::ValueType>(
                        (*i_chunk)->GetVariableData(),
                        std::distance(mData.begin(), i_chunk),
                        r_stored_history_policy.GetStepCategory());
                }
            }
        }

        KRATOS_ERROR << "No data found in container for " << rVariable.Name() << " with the requested policies" << std::endl;
    }

    /**
     * @brief Resize all sparse data chunks indexed by the given accessor.
     * @details Every sparse chunk whose policy references rIndexAccessor is rebuilt (its data is discarded and zero-initialized) with as many entities as there are active entries (index > -1) in the index span. Existing spans/accessors into the old chunk data are invalidated.
     * @param rIndexAccessor Accessor for a Variable<int> set to the sparse index for active entries or to a negative value otherwise.
     */
    void UpdateSparseStorage(const DataAccessor<int>& rIndexAccessor);

    /**
     * @brief Append entries for new entities to the sparse chunks indexed by the given accessor.
     * @details Entities of rEntityIndices whose index span entry is -1 get the next free sparse index assigned (updating the index span), and every sparse chunk referencing rIndexAccessor is grown accordingly, preserving its existing values and zero-initializing the appended entries (contrary to UpdateSparseStorage, which rebuilds from scratch). Spans obtained before the call are invalidated.
     * @param rIndexAccessor Accessor for the Variable<int> holding the sparse indices.
     * @param rEntityIndices The (entity) positions in the index span to add to the sparse storage.
     */
    void AddToSparseStorage(const DataAccessor<int>& rIndexAccessor, const std::vector<int>& rEntityIndices);

    /**
     * @brief Resize all dense (non-sparse) chunks to the given number of entities per step.
     * @details Delegates to DataChunkBase::ResizeData, which preserves the first min(old, new) values of every step slot and zero-initializes appended entries with the chunk's policy zero. Chunks already holding NewNumberOfEntities entities are left untouched, so repeated calls with the same size are cheap. Sparse chunks are skipped entirely: their sizing is governed by UpdateSparseStorage / AddToSparseStorage. All previously obtained spans over resized chunks are invalidated. Not thread-safe; external synchronization is required.
     * @param NewNumberOfEntities The new number of entities per step of the dense chunks.
     */
    void Resize(std::size_t NewNumberOfEntities);

    /**
     * @brief Clone the data of the current step onto the next step slot for all chunks of the given category.
     * @details Delegates to DataChunk::CloneStepData: only chunks whose history policy matches the category advance their ring and copy their data; all other chunks are untouched.
     * @param rStepCategory The step category being cloned.
     */
    void CloneStepData(StepCategory rStepCategory);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the data span of the current step for the given variable and policy.
     * @details The chunk is searched by variable and value policy (IsSame, i.e. the configuration must match too); the span covers the latest step slot. Raises an error if no matching chunk exists.
     * @tparam TVariableType The variable type.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable of the requested data.
     * @param rValuePolicy The value policy the data was added with.
     * @return A span over the values of the current step (one value per entity).
     */
    template<typename TVariableType, typename TValuePolicyType>
    typename TValuePolicyType::StepSpanType GetDataSpan(const TVariableType& rVariable, const TValuePolicyType& rValuePolicy)
    {
        for (const auto& p_chunk : mData) {
            if (p_chunk->GetVariableData() == rVariable && p_chunk->GetValuePolicy().IsSame(rValuePolicy)) {
                auto p_typed_chunk = std::dynamic_pointer_cast<DataChunk<typename TValuePolicyType::ValueType>>(p_chunk);
                if (p_typed_chunk) {
                    auto& r_history_policy = p_typed_chunk->GetHistoryPolicy();
                    const std::size_t current_step_index = r_history_policy.GetLatestStepIndex();
                    return typename TValuePolicyType::StepSpanType(p_typed_chunk->GetDataPointer(current_step_index),
                                                                   p_typed_chunk->NumberOfEntitiesPerStep());
                } else {
                    KRATOS_ERROR << "Failed to cast DataChunkBase to DataChunk<TValueType>." << std::endl;
                }
            }
        }

        // If the variable and policy are not found, throw an error
        KRATOS_ERROR << "Variable or policy not found in DataContainer for " << rVariable.Name() << "." << std::endl;
    }

    /**
     * @brief Get the data span of a specific step for the given variable and policy.
     * @details Like GetDataSpan(rVariable, rValuePolicy) but for the step slot given by the category and the number of steps before the current one (resolved by the chunk's history policy).
     * @tparam TVariableType The variable type.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable of the requested data.
     * @param rValuePolicy The value policy the data was added with.
     * @param rStepCategory The step category of the request.
     * @param rStepBeforeCurrent How many steps before the current one is requested.
     * @return A span over the values of the requested step (one value per entity).
     */
    template<typename TVariableType, typename TValuePolicyType>
    typename TValuePolicyType::StepSpanType GetDataSpan(const TVariableType& rVariable, const TValuePolicyType& rValuePolicy, StepCategory rStepCategory, int rStepBeforeCurrent = 0)
    {
        for (const auto& p_chunk : mData) {
            if (p_chunk->GetVariableData() == rVariable && p_chunk->GetValuePolicy().IsSame(rValuePolicy)) {
                auto p_typed_chunk = std::dynamic_pointer_cast<DataChunk<typename TValuePolicyType::ValueType>>(p_chunk);
                if (p_typed_chunk) {
                    auto& r_history_policy = p_typed_chunk->GetHistoryPolicy();
                    const std::size_t current_step_index = r_history_policy.GetStepIndex(rStepCategory, rStepBeforeCurrent);
                    return typename TValuePolicyType::StepSpanType(p_typed_chunk->GetDataPointer(current_step_index),
                                                                   p_typed_chunk->NumberOfEntitiesPerStep());
                } else {
                    KRATOS_ERROR << "Failed to cast DataChunkBase to DataChunk<TValueType>." << std::endl;
                }
            }
        }

        // If the variable and policy are not found, throw an error
        KRATOS_ERROR << "Variable or policy not found in DataContainer for " << rVariable.Name() << "." << std::endl;
    }

    /**
     * @brief Get the data span for the given accessor.
     * @details O(1): the chunk is addressed directly by the index cached in the accessor; the step slot is resolved from the accessor's category and step offset. The accessor must have been produced by this container (the variable identity is verified; an accessor from a different container raises an error).
     * @tparam TAccessorType The accessor type.
     * @param rAccessor The accessor of the requested data.
     * @return A span over the values of the accessor's step (one value per entity).
     */
    template<typename TAccessorType>
    typename TAccessorType::StepSpanType GetDataSpan(const TAccessorType& rAccessor)
    {
        KRATOS_ERROR_IF_NOT(rAccessor.GetIndex() < mData.size())
            << "Accessor index out of bounds: " << rAccessor.GetIndex() << ". DataContainer size: " << mData.size() << std::endl;

        auto& p_chunk = mData[rAccessor.GetIndex()];
        // check has the same variable. (This can happen when the accessor is from a different container)
        KRATOS_ERROR_IF_NOT(p_chunk->GetVariableData() == rAccessor.GetVariableData())
            << "Variable mismatch in DataContainer." << std::endl;

        // we should check if this affects the performance, if so we can use a static_cast and dynamic_cast only in full debug mode
        auto p_typed_chunk = std::dynamic_pointer_cast<DataChunk<typename TAccessorType::ValueType>>(p_chunk);
        if (p_typed_chunk) {
            auto& r_history_policy = p_typed_chunk->GetHistoryPolicy();
            const std::size_t current_step_index = r_history_policy.GetStepIndex(rAccessor.GetStepCategory(), rAccessor.GetStepBeforeCurrent());
            return typename TAccessorType::StepSpanType(p_typed_chunk->GetDataPointer(current_step_index),
                                                        p_typed_chunk->NumberOfEntitiesPerStep());
        } else {
            if (rAccessor.pGetVariableData() == nullptr)
                KRATOS_ERROR << "Accessor variable is null. The accessor is not initialized properly. Please check if the variable is added to this container" << std::endl;
            else if (rAccessor.GetVariableData().Name() != p_chunk->GetVariableData().Name())
                KRATOS_ERROR << "Accessor variable mismatch. Accessor variable: " << rAccessor.GetVariableData().Name()
                    << ", Chunk variable: " << p_chunk->GetVariableData().Name() << std::endl;
            else
                KRATOS_ERROR << "Failed to find the corresponding data for the given accessor. Please check if the accessor is for this container."
                    << "Accessor variable: " << rAccessor.GetVariableData().Name() << ", Chunk variable: " << p_chunk->GetVariableData().Name() << std::endl;
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the variable exists in the container.
     * @tparam TVariableType The variable type.
     * @param rVariable The variable to look for.
     * @return true if a chunk for the variable exists.
     */
    template<typename TVariableType>
    bool Has(const TVariableType& rVariable) const
    {
        for (const auto& p_chunk : mData) {
            if (p_chunk->GetVariableData() == rVariable) {
                return true; // Variable found
            }
        }
        return false; // Variable not found
    }

    /**
     * @brief Check if the variable exists in the container with the given value policy type.
     * @details Only the policy type is compared (IsSameType), not its configuration.
     * @tparam TVariableType The variable type.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable to look for.
     * @param rValuePolicy A value policy of the type to look for.
     * @return true if a chunk for the variable with the given value policy type exists.
     */
    template<typename TVariableType, typename TValuePolicyType>
    bool Has(const TVariableType& rVariable, const TValuePolicyType& rValuePolicy) const
    {
        for (const auto& p_chunk : mData) {
            if (p_chunk->GetVariableData() == rVariable && p_chunk->GetValuePolicy().IsSameType(rValuePolicy)) {
                return true; // Variable with the specific policy found
            }
        }
        return false; // Variable with the specific policy not found
    }

    /**
     * @brief Check if the variable exists in the container with the given value and history policy types.
     * @details Only the policy types are compared (IsSameType), not their configuration.
     * @tparam TVariableType The variable type.
     * @param rVariable The variable to look for.
     * @param rValuePolicy A value policy of the type to look for.
     * @param rHistoryPolicy A history policy of the type to look for.
     * @return true if a chunk for the variable with the given policy types exists.
     */
    template<typename TVariableType>
    bool Has(const TVariableType& rVariable, const DataValuePolicyBase& rValuePolicy, const DataHistoryPolicyBase& rHistoryPolicy) const
    {
        for (const auto& p_chunk : mData) {
            if (p_chunk->GetVariableData() == rVariable &&
                p_chunk->GetValuePolicy().IsSameType(rValuePolicy) &&
                p_chunk->GetHistoryPolicy().IsSameType(rHistoryPolicy)) {
                return true; // Variable with the specific policies found
            }
        }
        return false; // Variable with the specific policies not found
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<std::shared_ptr<DataChunkBase>> mData; /// The data stored in the container. Each element is a shared pointer to the chunk

    std::size_t mChunkSize = 256;                      /// Number of entities in each dense chunk

    LockObject mMutex;                                 /// Lock for mutual exclusion when adding data chunks to the container

    ///@}

}; // Class DataContainer

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DataContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos
