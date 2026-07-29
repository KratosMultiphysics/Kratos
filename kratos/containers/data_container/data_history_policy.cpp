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
#include "includes/define.h"
#include "containers/data_container/data_history_policy.h"

namespace Kratos
{

std::unique_ptr<DataHistoryPolicyBase> NonHistoricalDataPolicy::Clone() const
{
    return std::make_unique<NonHistoricalDataPolicy>();
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t NonHistoricalDataPolicy::CloneStepIndex(StepCategory Category)
{
    return 0; // Non-historical data policy does not support cloning steps
}

/***********************************************************************************/
/***********************************************************************************/

bool NonHistoricalDataPolicy::IsSameType(const DataHistoryPolicyBase& rOther) const
{
    return dynamic_cast<const NonHistoricalDataPolicy*>(&rOther) != nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

bool NonHistoricalDataPolicy::IsSame(const DataHistoryPolicyBase& rOther) const
{
    return IsSameType(rOther);
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t NonHistoricalDataPolicy::GetLatestStepIndex() const
{
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t NonHistoricalDataPolicy::GetStepIndex(StepCategory Category, int /*StepBeforeCurrent*/) const
{
    KRATOS_ERROR_IF(Category != StepCategory::AnyStep) << "NonHistoricalDataPolicy does not support step categories." << std::endl;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

StepCategory NonHistoricalDataPolicy::GetStepCategory() const
{
    return StepCategory::AnyStep;
}

/***********************************************************************************/
/***********************************************************************************/

HistoricalDataPolicy::HistoricalDataPolicy(StepCategory TheCategory, std::size_t TotalNumberOfSteps)
    : DataHistoryPolicyBase(TotalNumberOfSteps),
      mCategory(TheCategory),
      mLatestStepIndex(0)
{
    KRATOS_DEBUG_ERROR_IF(TotalNumberOfSteps == 0) << "Total number of steps must be greater than zero." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::unique_ptr<DataHistoryPolicyBase> HistoricalDataPolicy::Clone() const
{
    return std::make_unique<HistoricalDataPolicy>(mCategory, mTotalNumberOfSteps);
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t HistoricalDataPolicy::CloneStepIndex(StepCategory Category)
{
    if (mCategory != Category) {
        return mLatestStepIndex; // No cloning if the category does not match
    }

    const std::size_t next_step_index = (mLatestStepIndex + 1) % mTotalNumberOfSteps;
    mLatestStepIndex = next_step_index;
    return next_step_index;
}

/***********************************************************************************/
/***********************************************************************************/

bool HistoricalDataPolicy::IsSameType(const DataHistoryPolicyBase& rOther) const
{
    return dynamic_cast<const HistoricalDataPolicy*>(&rOther) != nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

bool HistoricalDataPolicy::IsSame(const DataHistoryPolicyBase& rOther) const
{
    const auto* p_typed_other = dynamic_cast<const HistoricalDataPolicy*>(&rOther);
    if (p_typed_other) {
        return mTotalNumberOfSteps == p_typed_other->mTotalNumberOfSteps && mCategory == p_typed_other->mCategory;
    }
    return false; // Different type
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t HistoricalDataPolicy::GetLatestStepIndex() const
{
    return mLatestStepIndex;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t HistoricalDataPolicy::GetStepIndex(StepCategory Category, int StepBeforeCurrent) const
{
    KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(Category) >= static_cast<std::size_t>(StepCategory::NumberOfCategories))
        << "Invalid step category." << std::endl;
    KRATOS_DEBUG_ERROR_IF(StepBeforeCurrent < 0)
        << "StepBeforeCurrent must be non-negative." << std::endl;

    const std::size_t steps_before_current = static_cast<std::size_t>(StepBeforeCurrent);
    KRATOS_DEBUG_ERROR_IF(steps_before_current >= mTotalNumberOfSteps)
        << "StepBeforeCurrent must be less than the total number of steps." << std::endl;

    const std::size_t index = (mLatestStepIndex >= steps_before_current) ?
                              (mLatestStepIndex - steps_before_current) :
                              (mTotalNumberOfSteps + mLatestStepIndex - steps_before_current);

    return index;
}

/***********************************************************************************/
/***********************************************************************************/

StepCategory HistoricalDataPolicy::GetStepCategory() const
{
    return mCategory;
}

} // namespace Kratos
