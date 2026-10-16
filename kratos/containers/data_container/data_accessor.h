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
#include <span/span.hpp>

// Project includes
#include "containers/variable_data.h"
#include "containers/data_container/step_category.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DataAccessor
 * @ingroup KratosCore
 * @brief Lightweight, non-owning handle to the data of one variable stored in a DataContainer.
 * @details The accessor caches the index of the corresponding DataChunk inside its DataContainer, so the data can be retrieved in O(1) via DataContainer::GetDataSpan without searching the container by variable. It additionally stores which step of the chunk history it refers to (a StepCategory plus a number of steps before the current one), so re-parameterized accessors for other steps can be derived cheaply with GetStepAccessor.
 *
 * Accessors are produced by DataContainer::Add and DataContainer::GetAccessor; they hold a raw pointer to the canonical (Registry-owned) variable and do not own any data or keep the container alive. An accessor is only valid for the container that created it and as long as that container exists; DataContainer::GetDataSpan verifies the variable identity before returning data.
 *
 * @note This class is trivially copyable and cheap to pass by value. It provides no thread-safety guarantees of its own (it is immutable after construction).
 * @tparam TValueType Type of the values stored in the referenced DataChunk.
 * @author Pooyan Dadvand
 */
template<typename TValueType>
class DataAccessor
{
public:
    ///@name Type Definitions
    ///@{

    /// Type of the data stored in the chunk
    using ValueType = TValueType;

    /// Type of the span of data for a specific step
    using StepSpanType = Kratos::span<TValueType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Creates an uninitialized accessor (null variable, index 0, StepCategory::AnyStep). Using it with DataContainer::GetDataSpan raises an error.
     */
    DataAccessor()
        : mpVariable(nullptr),
          mIndex(0),
          mStepCategory(StepCategory::AnyStep),
          mStepBeforeCurrent(0)
    {}

    /**
     * @brief Constructor.
     * @param rVariable The variable associated with the data chunk. The accessor stores a pointer to it, so it must outlive the accessor (DataContainer passes the Registry-owned canonical variable, which lives for the whole program).
     * @param Index The index of the data chunk in the container.
     * @param StepCat The step category the accessor refers to.
     * @param StepBeforeCurrent How many steps before the current one the accessor refers to.
     */
    template<typename TVariableType>
    DataAccessor(
        TVariableType const& rVariable,
        std::size_t Index,
        StepCategory StepCat = StepCategory::TimeStep,
        int StepBeforeCurrent = 0)
        : mpVariable(&rVariable),
          mIndex(Index),
          mStepCategory(StepCat),
          mStepBeforeCurrent(StepBeforeCurrent)
    {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Equality comparison.
     * @details Two accessors are equal if they point to the same variable object (pointer identity, which holds for accessors produced by the same container thanks to the Registry canonicalization in DataContainer::Add), reference the same chunk index and the same step (category and steps before current).
     * @param rOther The accessor to compare against.
     * @return true if both accessors refer to the same data and step.
     */
    bool operator==(const DataAccessor<TValueType>& rOther) const
    {
        return (
            mpVariable == rOther.mpVariable &&
            mIndex == rOther.mIndex &&
            mStepCategory == rOther.mStepCategory &&
            mStepBeforeCurrent == rOther.mStepBeforeCurrent
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get an accessor to the same data but a different step.
     * @param StepCat The step category of the new accessor.
     * @param StepBeforeCurrent How many steps before the current one the new accessor refers to.
     * @return A new accessor for the same variable and chunk index, re-parameterized to the given step.
     */
    DataAccessor<TValueType> GetStepAccessor(StepCategory StepCat, int StepBeforeCurrent = 0) const
    {
        return DataAccessor<TValueType>(*mpVariable, mIndex, StepCat, StepBeforeCurrent);
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the index of the data chunk in the container.
     * @return The chunk index cached by this accessor.
     */
    std::size_t GetIndex() const
    {
        return mIndex;
    }

    /**
     * @brief Get the step category this accessor refers to.
     * @return The step category.
     */
    StepCategory GetStepCategory() const
    {
        return mStepCategory;
    }

    /**
     * @brief Get how many steps before the current one this accessor refers to.
     * @return The number of steps before the current one.
     */
    int GetStepBeforeCurrent() const
    {
        return mStepBeforeCurrent;
    }

    /**
     * @brief Get the variable associated with the data chunk.
     * @return Reference to the variable data. Must not be called on a default-constructed accessor.
     */
    const VariableData& GetVariableData() const
    {
        return *mpVariable;
    }

    /**
     * @brief Get the pointer to the variable associated with the data chunk.
     * @return Pointer to the variable data (nullptr for a default-constructed accessor).
     */
    const VariableData* pGetVariableData() const
    {
        return mpVariable;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    const VariableData* mpVariable;  /// Pointer to the variable associated with the data chunk

    std::size_t mIndex;              /// The index of the data chunk in the container

    StepCategory mStepCategory = StepCategory::AnyStep;  /// The step category for the data chunk (AnyStep when default-constructed, TimeStep by default otherwise)

    int mStepBeforeCurrent = 0;      /// The step before current for the data chunk, default is 0

    ///@}

}; // Class DataAccessor

///@}

///@} addtogroup block

} // namespace Kratos
