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
#include <concepts>
#include <unordered_map>

// External includes

// Project includes
#include "containers/data_container/data_container.h"

namespace Kratos::Future
{

///@name Type Definitions
///@{

/// Concept satisfied by Kratos entities providing an Id (Node, Element, Condition, Geometry, MasterSlaveConstraint, ...); keeps the entity-taking overloads from hijacking plain integral entity ids
template<class TEntityType>
concept EntityWithId = requires(const TEntityType& rEntity) {
    { rEntity.Id() } -> std::convertible_to<std::size_t>;
};

///@}
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class EntityDataContainer
 * @ingroup KratosCore
 * @brief Experimental side-car storage mapping Kratos entities onto a DataContainer.
 * @details This class is the Phase II bridge between Id-addressed Kratos entities (Node, Element, Condition, Geometry, MasterSlaveConstraint) and the slot-addressed, struct-of-arrays DataContainer: it owns one DataContainer plus a map from entity Id to a dense slot index. Entities are registered (RegisterEntity / RegisterEntityId, assigning slots in registration order) and variables are added per container (AddVariable for non-historical data, AddHistoricalVariable for TimeStep-buffered data, or the raw Add passthrough for sparse/custom policies); a value is then addressed as (variable chunk, entity slot).
 *
 * This is a PARALLEL storage path: it never touches, and is never touched by, the entities' existing storage (NodalData / DataValueContainer). Entities are identified purely by Id() — the entity classes are not modified and hold no reference to this container.
 *
 * Semantics and limitations (Phase II):
 * - Registration order defines slots; unregistering an entity only removes its Id from the map, leaving a hole (the slot is not reclaimed and the stored values remain in memory until the container is rebuilt).
 * - Dense chunks grow automatically (DataContainer::Resize) as entities register; growth invalidates previously obtained spans (and NumPy views).
 * - The buffer size for historical variables is fixed at construction.
 * - Only Registry-registered variables can be added (Phase I requirement).
 * - Copies are shallow: they share the chunk storage (DataContainer semantics) but copy the Id map.
 * - Thread-safety: registration and variable addition require external synchronization; afterwards, concurrent reads/writes to distinct entity slots are safe (distinct memory).
 * - No serialization and no MPI synchronization yet.
 * @see DataContainer, DataAccessor, ModelPartDataContainer
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) EntityDataContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EntityDataContainer
    KRATOS_CLASS_POINTER_DEFINITION(EntityDataContainer);

    /// Type of the entity identifiers (entity Id)
    using IndexType = std::size_t;

    /// Type of the dense slot indices inside the chunks
    using SlotType = std::size_t;

    /// Accessor type alias for the stored data of a given value type
    template<typename TValueType>
    using DataAccessor = Kratos::DataAccessor<TValueType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param BufferSize Number of steps stored by historical variables (including the current one). Fixed for the lifetime of the container.
     * @param ChunkSize Initial number of entity slots of the dense chunks.
     */
    explicit EntityDataContainer(std::size_t BufferSize = 1, std::size_t ChunkSize = 256);

    /// Destructor.
    ~EntityDataContainer() = default;

    /**
     * @brief Copy constructor.
     * @details Shallow with respect to the data (chunks are shared with rOther, DataContainer semantics); the Id map and variable records are copied.
     * @param rOther The container to copy from.
     */
    EntityDataContainer(const EntityDataContainer& rOther) = default;

    /// Move constructor.
    EntityDataContainer(EntityDataContainer&& rOther) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy assignment operator (shallow data sharing, see copy constructor).
    EntityDataContainer& operator=(const EntityDataContainer& rOther) = default;

    /// Move assignment operator.
    EntityDataContainer& operator=(EntityDataContainer&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Register an entity Id, assigning it a dense slot.
     * @details Idempotent: registering an already known Id returns its existing slot. The dense chunks grow automatically when the number of slots exceeds the current capacity (invalidating previously obtained spans).
     * @param EntityId The Id of the entity to register.
     * @return The slot assigned to the entity.
     */
    SlotType RegisterEntityId(IndexType EntityId);

    /**
     * @brief Register an entity, assigning it a dense slot.
     * @details Convenience overload of RegisterEntityId using rEntity.Id().
     * @tparam TEntityType Any type providing an Id() method (Node, Element, Condition, Geometry, MasterSlaveConstraint, ...).
     * @param rEntity The entity to register.
     * @return The slot assigned to the entity.
     */
    template<EntityWithId TEntityType>
    SlotType RegisterEntity(const TEntityType& rEntity)
    {
        return RegisterEntityId(rEntity.Id());
    }

    /**
     * @brief Register a range of entities in iteration order with a single storage growth.
     * @tparam TIteratorType Iterator dereferencing to an entity providing Id().
     * @param First Begin of the range.
     * @param Last End of the range.
     */
    template<class TIteratorType>
    void RegisterEntities(TIteratorType First, TIteratorType Last)
    {
        KRATOS_TRY

        for (auto it = First; it != Last; ++it) {
            const IndexType entity_id = it->Id();
            if (mIdToSlot.find(entity_id) == mIdToSlot.end()) {
                mIdToSlot.emplace(entity_id, mNextSlot++);
            }
        }
        EnsureCapacity();

        KRATOS_CATCH("")
    }

    /**
     * @brief Unregister an entity Id.
     * @details Phase II limitation: only the Id-to-slot entry is removed. The slot becomes a hole — it is not reclaimed, the span length does not shrink, and the stored values remain until the container is rebuilt. Re-registering the same Id assigns a NEW slot.
     * @param EntityId The Id of the entity to unregister.
     */
    void UnregisterEntityId(IndexType EntityId);

    /**
     * @brief Add a non-historical dense variable (one value per entity).
     * @details Maps onto DataContainer::Add with DataValuePolicy<T>(rVariable.Zero()) and NonHistoricalDataPolicy. The variable must be registered in the Kratos Registry.
     * @tparam TValueType The value type of the variable.
     * @param rVariable The variable to add.
     * @return An accessor to the data (index it with Index(entity) inside GetDataSpan).
     */
    template<typename TValueType>
    DataAccessor<TValueType> AddVariable(const Variable<TValueType>& rVariable)
    {
        KRATOS_TRY

        return Add(rVariable, DataValuePolicy<TValueType>(rVariable.Zero()), NonHistoricalDataPolicy());

        KRATOS_CATCH("")
    }

    /**
     * @brief Add a historical dense variable buffered over the container's buffer size.
     * @details Maps onto DataContainer::Add with DataValuePolicy<T>(rVariable.Zero()) and HistoricalDataPolicy(StepCategory::TimeStep, GetBufferSize()) — the DataContainer equivalent of a nodal solution-step variable. Advance the buffer with CloneStepData(StepCategory::TimeStep).
     * @tparam TValueType The value type of the variable.
     * @param rVariable The variable to add.
     * @return An accessor to the current-step data.
     */
    template<typename TValueType>
    DataAccessor<TValueType> AddHistoricalVariable(const Variable<TValueType>& rVariable)
    {
        KRATOS_TRY

        return Add(rVariable, DataValuePolicy<TValueType>(rVariable.Zero()), HistoricalDataPolicy(StepCategory::TimeStep, mBufferSize));

        KRATOS_CATCH("")
    }

    /**
     * @brief Add a variable with full control over the value and history policies.
     * @details Raw passthrough to DataContainer::Add (e.g. for SparseDataValuePolicy or custom zeros/categories). Phase I semantics apply, including the caller's responsibility for a meaningful policy zero.
     * @tparam TVariableValueType The value type of the variable.
     * @tparam TValuePolicyType The concrete value policy type.
     * @param rVariable The variable to add.
     * @param rValuePolicy The value policy (cloned into the chunk).
     * @param rHistoryPolicy The history policy (cloned into the chunk).
     * @return An accessor to the data.
     */
    template<typename TVariableValueType, typename TValuePolicyType>
    DataAccessor<typename TValuePolicyType::ValueType> Add(
        const Variable<TVariableValueType>& rVariable,
        const TValuePolicyType& rValuePolicy,
        const DataHistoryPolicyBase& rHistoryPolicy = NonHistoricalDataPolicy())
    {
        KRATOS_TRY

        auto accessor = mDataContainer.Add(rVariable, rValuePolicy, rHistoryPolicy);
        mVariableRecords[rVariable.Key()] = VariableRecord{accessor.GetIndex(), rHistoryPolicy.GetStepCategory()};
        // A chunk created after earlier growth starts at the construction chunk size; bring it to capacity
        mDataContainer.Resize(mCapacity);
        return accessor;

        KRATOS_CATCH("")
    }

    /**
     * @brief Clone the current step onto the next step slot for all chunks of the given category.
     * @details Passthrough to DataContainer::CloneStepData — the DataContainer equivalent of Node::CloneSolutionStepData; call it right after ModelPart::CloneTimeStep when mirroring the legacy workflow.
     * @param Category The step category being cloned.
     */
    void CloneStepData(StepCategory Category);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the slot of a registered entity Id.
     * @param EntityId The Id of the entity.
     * @return The dense slot of the entity. Raises an error for unknown Ids.
     */
    SlotType Index(IndexType EntityId) const;

    /**
     * @brief Get the slot of a registered entity.
     * @tparam TEntityType Any type providing an Id() method.
     * @param rEntity The entity.
     * @return The dense slot of the entity. Raises an error for unregistered entities.
     */
    template<EntityWithId TEntityType>
    SlotType Index(const TEntityType& rEntity) const
    {
        return Index(rEntity.Id());
    }

    /**
     * @brief Get a reference to the value of a variable for one entity.
     * @details O(1): rebuilds the accessor from the record stored at Add time and indexes the step span at the entity slot. For historical variables StepBeforeCurrent selects the step (0 = current); for non-historical variables it must be 0.
     * @tparam TValueType The value type of the variable.
     * @param EntityId The Id of the entity.
     * @param rVariable The variable (must have been added).
     * @param StepBeforeCurrent How many steps before the current one (historical variables only).
     * @return Reference to the stored value.
     */
    template<typename TValueType>
    TValueType& GetValue(IndexType EntityId, const Variable<TValueType>& rVariable, int StepBeforeCurrent = 0)
    {
        KRATOS_TRY

        const auto& r_record = GetVariableRecord(rVariable, StepBeforeCurrent);
        DataAccessor<TValueType> accessor(rVariable, r_record.ChunkIndex, r_record.Category, StepBeforeCurrent);
        return mDataContainer.GetDataSpan(accessor)[Index(EntityId)];

        KRATOS_CATCH("")
    }

    /**
     * @brief Get a reference to the value of a variable for one entity.
     * @tparam TEntityType Any type providing an Id() method.
     * @tparam TValueType The value type of the variable.
     * @param rEntity The entity.
     * @param rVariable The variable (must have been added).
     * @param StepBeforeCurrent How many steps before the current one (historical variables only).
     * @return Reference to the stored value.
     */
    template<EntityWithId TEntityType, typename TValueType>
    TValueType& GetValue(const TEntityType& rEntity, const Variable<TValueType>& rVariable, int StepBeforeCurrent = 0)
    {
        return GetValue(rEntity.Id(), rVariable, StepBeforeCurrent);
    }

    /**
     * @brief Set the value of a variable for one entity.
     * @tparam TValueType The value type of the variable.
     * @param EntityId The Id of the entity.
     * @param rVariable The variable (must have been added).
     * @param rValue The value to assign.
     * @param StepBeforeCurrent How many steps before the current one (historical variables only).
     */
    template<typename TValueType>
    void SetValue(IndexType EntityId, const Variable<TValueType>& rVariable, const TValueType& rValue, int StepBeforeCurrent = 0)
    {
        GetValue(EntityId, rVariable, StepBeforeCurrent) = rValue;
    }

    /**
     * @brief Set the value of a variable for one entity.
     * @tparam TEntityType Any type providing an Id() method.
     * @tparam TValueType The value type of the variable.
     * @param rEntity The entity.
     * @param rVariable The variable (must have been added).
     * @param rValue The value to assign.
     * @param StepBeforeCurrent How many steps before the current one (historical variables only).
     */
    template<EntityWithId TEntityType, typename TValueType>
    void SetValue(const TEntityType& rEntity, const Variable<TValueType>& rVariable, const TValueType& rValue, int StepBeforeCurrent = 0)
    {
        GetValue(rEntity.Id(), rVariable, StepBeforeCurrent) = rValue;
    }

    /**
     * @brief Get an accessor for a previously added variable.
     * @details The fast bulk-access pattern is GetDataSpan(accessor)[Index(entity)] inside loops.
     * @tparam TValueType The value type of the variable.
     * @param rVariable The variable (must have been added).
     * @return An accessor addressing the variable's chunk with the category recorded at Add time.
     */
    template<typename TValueType>
    DataAccessor<TValueType> GetAccessor(const Variable<TValueType>& rVariable) const
    {
        KRATOS_TRY

        const auto& r_record = GetVariableRecord(rVariable, 0);
        return DataAccessor<TValueType>(rVariable, r_record.ChunkIndex, r_record.Category);

        KRATOS_CATCH("")
    }

    /**
     * @brief Get the data span for the given accessor.
     * @details Passthrough to DataContainer::GetDataSpan: one value per slot (index with Index(entity)). Spans are invalidated by entity registration growth and Resize.
     * @tparam TAccessorType The accessor type.
     * @param rAccessor The accessor of the requested data.
     * @return A span over the values of the accessor's step.
     */
    template<typename TAccessorType>
    typename TAccessorType::StepSpanType GetDataSpan(const TAccessorType& rAccessor)
    {
        return mDataContainer.GetDataSpan(rAccessor);
    }

    /**
     * @brief Get the underlying DataContainer.
     * @details Escape hatch for advanced use (e.g. sparse storage updates, NumPy span access through the Phase I bindings). Slot bookkeeping remains this class's responsibility.
     * @return Reference to the owned DataContainer.
     */
    DataContainer& GetDataContainer();

    /**
     * @brief Get the number of steps stored by historical variables.
     * @return The buffer size fixed at construction.
     */
    std::size_t GetBufferSize() const;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if an entity Id is registered.
     * @param EntityId The Id of the entity.
     * @return true if the Id has a slot.
     */
    bool HasEntity(IndexType EntityId) const;

    /**
     * @brief Check if a variable has been added to this container.
     * @param rVariable The variable to look for.
     * @return true if the variable was added.
     */
    bool Has(const VariableData& rVariable) const;

    /**
     * @brief Get the number of registered entities.
     * @return The number of currently registered Ids (holes from unregistration do not count).
     */
    std::size_t NumberOfEntities() const;

    /**
     * @brief Get the number of slots per step of the dense chunks.
     * @details This is the span length; it only grows (holes from unregistration are not reclaimed).
     * @return The current slot capacity.
     */
    std::size_t Capacity() const;

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
    ///@name Type Definitions
    ///@{

    /// Bookkeeping of one added variable: which chunk it lives in and its step category
    struct VariableRecord
    {
        std::size_t ChunkIndex;
        StepCategory Category;
    };

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Grow the dense chunks so that every assigned slot fits.
     * @details Doubles the capacity (at least to the highest assigned slot count) and resizes the dense chunks, preserving existing values. No-op while the slots fit.
     */
    void EnsureCapacity();

    /**
     * @brief Get the record of an added variable, validating the step request.
     * @param rVariable The variable to look up.
     * @param StepBeforeCurrent The requested step offset (validated against the record's category).
     * @return The variable record. Raises clear errors for unknown variables and for step requests on non-historical variables.
     */
    const VariableRecord& GetVariableRecord(const VariableData& rVariable, int StepBeforeCurrent) const;

    ///@}
    ///@name Member Variables
    ///@{

    DataContainer mDataContainer;                                        /// The chunked storage (one chunk per added variable)

    std::unordered_map<IndexType, SlotType> mIdToSlot;                   /// Map from entity Id to dense slot

    std::unordered_map<VariableData::KeyType, VariableRecord> mVariableRecords; /// Per-variable chunk index and step category, keyed by Variable::Key()

    SlotType mNextSlot = 0;                                              /// The next slot to assign

    std::size_t mCapacity;                                               /// Current number of slots per step of the dense chunks

    std::size_t mBufferSize;                                             /// Number of steps stored by historical variables

    ///@}

}; // Class EntityDataContainer

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const EntityDataContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos::Future
