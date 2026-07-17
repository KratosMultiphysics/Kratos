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
#include <cstddef>
#include <memory>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_container/step_category.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DataHistoryPolicyBase
 * @ingroup KratosCore
 * @brief Interface defining how the step history of the data is stored in each DataChunk.
 * @details A history policy determines how many steps of data a DataChunk keeps (GetTotalNumberOfSteps) and maps step requests (a StepCategory plus a number of steps before the current one) to the index of the corresponding step slot inside the chunk storage. Historical policies implement a ring buffer over the step slots: cloning a step (CloneStepIndex) advances the latest-step index only when the requested category matches the policy's own category, so containers holding chunks of mixed categories can be cloned selectively per category.
 *
 * Policies are cloned into each DataChunk; HistoricalDataPolicy is stateful (it tracks the latest step index) and provides no thread-safety guarantees.
 * @see NonHistoricalDataPolicy, HistoricalDataPolicy, DataChunk, DataContainer
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) DataHistoryPolicyBase
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param TotalNumberOfSteps Number of steps to be stored including the current one.
     */
    explicit DataHistoryPolicyBase(std::size_t TotalNumberOfSteps = 1)
        : mTotalNumberOfSteps(TotalNumberOfSteps) {}

    /// Destructor.
    virtual ~DataHistoryPolicyBase() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the history policy.
     * @return A newly allocated copy of this policy. Note stateful policies reset their running step index (see HistoricalDataPolicy::Clone).
     */
    virtual std::unique_ptr<DataHistoryPolicyBase> Clone() const = 0;

    /**
     * @brief Advance the step ring buffer for the given category.
     * @details Only advances (and returns the new latest index) when the given category matches the policy's own category; otherwise the latest step index is returned unchanged. This is what makes DataContainer::CloneStepData act selectively on the chunks of the requested category.
     * @param Category The step category being cloned.
     * @return The step index that became (or remains) the latest one.
     */
    virtual std::size_t CloneStepIndex(StepCategory Category) = 0;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the given policy has the same dynamic type as this one.
     * @details Configuration details (category, number of steps) are ignored; see IsSame.
     * @param rOther The policy to compare against.
     * @return true if both policies are of the same dynamic type.
     */
    virtual bool IsSameType(const DataHistoryPolicyBase& rOther) const = 0;

    /**
     * @brief Check if the given policy is the same as this one (type and configuration).
     * @param rOther The policy to compare against.
     * @return true if both policies are of the same dynamic type with identical configuration.
     */
    virtual bool IsSame(const DataHistoryPolicyBase& rOther) const = 0;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the total number of steps to be stored (including the current one).
     * @return The total number of steps.
     */
    std::size_t GetTotalNumberOfSteps() const
    {
        return mTotalNumberOfSteps;
    }

    /**
     * @brief Set the total number of steps to be stored (including the current one).
     * @param TotalNumberOfSteps The new total number of steps.
     */
    void SetTotalNumberOfSteps(std::size_t TotalNumberOfSteps)
    {
        mTotalNumberOfSteps = TotalNumberOfSteps;
    }

    /**
     * @brief Get the index of the latest (current) step slot.
     * @return The latest step index.
     */
    virtual std::size_t GetLatestStepIndex() const = 0;

    /**
     * @brief Get the step slot index for a specific category and step offset.
     * @param Category The step category of the request.
     * @param StepBeforeCurrent How many steps before the current one is requested.
     * @return The index of the corresponding step slot.
     */
    virtual std::size_t GetStepIndex(StepCategory Category, int StepBeforeCurrent = 0) const = 0;

    /**
     * @brief Get the step category tracked by this policy.
     * @return The step category.
     */
    virtual StepCategory GetStepCategory() const = 0;

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    std::size_t mTotalNumberOfSteps; /// Number of steps to be stored including the current one

    ///@}

}; // Class DataHistoryPolicyBase

/**
 * @class NonHistoricalDataPolicy
 * @ingroup KratosCore
 * @brief History policy storing only one value (no history) per entity.
 * @details The single stored step always has index 0 and belongs to StepCategory::AnyStep. Requesting a step index for any other category is an error, and cloning never advances anything.
 * @see DataHistoryPolicyBase, HistoricalDataPolicy
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) NonHistoricalDataPolicy : public DataHistoryPolicyBase
{
public:
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NonHistoricalDataPolicy() : DataHistoryPolicyBase(1) {}

    /// Destructor.
    ~NonHistoricalDataPolicy() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the history policy.
     * @return A newly allocated NonHistoricalDataPolicy.
     */
    std::unique_ptr<DataHistoryPolicyBase> Clone() const override;

    /**
     * @brief Advance the step ring buffer for the given category.
     * @details Non-historical data has no history to clone; the single step index 0 is always returned.
     * @return Always 0.
     */
    std::size_t CloneStepIndex(StepCategory Category) override;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the given policy has the same dynamic type as this one.
     * @param rOther The policy to compare against.
     * @return true if rOther is a NonHistoricalDataPolicy.
     */
    bool IsSameType(const DataHistoryPolicyBase& rOther) const override;

    /**
     * @brief Check if the given policy is the same as this one.
     * @param rOther The policy to compare against.
     * @return true if rOther is a NonHistoricalDataPolicy (the policy has no configuration).
     */
    bool IsSame(const DataHistoryPolicyBase& rOther) const override;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the index of the latest (current) step slot.
     * @return Always 0 (non-historical data has only one step, the current one).
     */
    std::size_t GetLatestStepIndex() const override;

    /**
     * @brief Get the step slot index for a specific category and step offset.
     * @param Category The step category of the request. Must be StepCategory::AnyStep.
     * @return Always 0.
     */
    std::size_t GetStepIndex(StepCategory Category, int StepBeforeCurrent = 0) const override;

    /**
     * @brief Get the step category tracked by this policy.
     * @return StepCategory::AnyStep (non-historical data has no specific category).
     */
    StepCategory GetStepCategory() const override;

    ///@}

}; // Class NonHistoricalDataPolicy

/**
 * @class HistoricalDataPolicy
 * @ingroup KratosCore
 * @brief History policy storing data for each step, up to a maximum number of steps, in a ring buffer.
 * @details The policy is bound to one StepCategory (e.g. time steps, iteration steps, sub-steps) and keeps track of the latest step slot index. Cloning advances the ring (wrapping around at the total number of steps) only when the cloned category matches the policy's own category. GetStepIndex maps "N steps before the current one" onto the ring, wrapping backwards.
 * @note Cloning the policy object itself (Clone) resets the running latest-step index to zero: freshly created chunks (DataChunk::CreateNew, DataContainer::Initialize) start their history at slot 0.
 * @see DataHistoryPolicyBase, NonHistoricalDataPolicy
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) HistoricalDataPolicy : public DataHistoryPolicyBase
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param TheCategory Category of the step tracked by this policy (e.g. TimeStep, IterationStep).
     * @param TotalNumberOfSteps Number of steps to store including the current one. Must be greater than zero.
     */
    explicit HistoricalDataPolicy(StepCategory TheCategory, std::size_t TotalNumberOfSteps);

    /// Destructor.
    ~HistoricalDataPolicy() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the history policy.
     * @return A newly allocated HistoricalDataPolicy with the same category and number of steps, but with the latest step index reset to zero.
     */
    std::unique_ptr<DataHistoryPolicyBase> Clone() const override;

    /**
     * @brief Advance the step ring buffer for the given category.
     * @details Advances the latest step index (wrapping around at the total number of steps) only when the given category matches this policy's category; otherwise the index is left untouched.
     * @param Category The step category being cloned.
     * @return The step index that became (or remains) the latest one.
     */
    std::size_t CloneStepIndex(StepCategory Category) override;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the given policy has the same dynamic type as this one.
     * @param rOther The policy to compare against.
     * @return true if rOther is a HistoricalDataPolicy (category and number of steps are ignored).
     */
    bool IsSameType(const DataHistoryPolicyBase& rOther) const override;

    /**
     * @brief Check if the given policy is the same as this one.
     * @param rOther The policy to compare against.
     * @return true if rOther is a HistoricalDataPolicy with the same category and total number of steps.
     */
    bool IsSame(const DataHistoryPolicyBase& rOther) const override;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the index of the latest (current) step slot.
     * @return The latest step index.
     */
    std::size_t GetLatestStepIndex() const override;

    /**
     * @brief Get the step slot index for a specific category and step offset.
     * @details Maps "StepBeforeCurrent steps before the latest one" onto the ring buffer, wrapping backwards over the total number of steps.
     * @param Category The step category of the request.
     * @param StepBeforeCurrent How many steps before the current one is requested. Must be non-negative and less than the total number of steps.
     * @return The index of the corresponding step slot.
     */
    std::size_t GetStepIndex(StepCategory Category, int StepBeforeCurrent = 0) const override;

    /**
     * @brief Get the step category tracked by this policy.
     * @return The step category.
     */
    StepCategory GetStepCategory() const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    StepCategory mCategory;           /// Category of the step (e.g., TimeStep, IterationStep, etc.)

    std::size_t mLatestStepIndex = 0; /// Index of the latest step added

    ///@}

}; // Class HistoricalDataPolicy

///@}

///@} addtogroup block

} // namespace Kratos
