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
#include "containers/data_container/data_history_policy.h"
#include "includes/expect.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataHistoryPolicyComparison, KratosCoreFastSuite)
{
    NonHistoricalDataPolicy non_historical_policy;

    KRATOS_EXPECT_TRUE(NonHistoricalDataPolicy().IsSameType(non_historical_policy));
    KRATOS_EXPECT_TRUE(NonHistoricalDataPolicy().IsSame(non_historical_policy));

    HistoricalDataPolicy historical_policy(StepCategory::TimeStep, 3);

    KRATOS_EXPECT_FALSE(historical_policy.IsSameType(non_historical_policy));
    KRATOS_EXPECT_FALSE(non_historical_policy.IsSameType(historical_policy));

    HistoricalDataPolicy other_historical_policy(StepCategory::TimeStep, 2);

    KRATOS_EXPECT_TRUE(historical_policy.IsSameType(other_historical_policy)); // same type
    KRATOS_EXPECT_FALSE(historical_policy.IsSame(other_historical_policy)); // different number of steps

    HistoricalDataPolicy different_category_policy(StepCategory::IterationStep, 3);
    KRATOS_EXPECT_TRUE(historical_policy.IsSameType(different_category_policy)); // same type
    KRATOS_EXPECT_FALSE(historical_policy.IsSame(different_category_policy)); // different category

    HistoricalDataPolicy equal_policy(StepCategory::TimeStep, 3);
    KRATOS_EXPECT_TRUE(historical_policy.IsSame(equal_policy)); // same category and number of steps
}

KRATOS_TEST_CASE_IN_SUITE(NonHistoricalDataPolicyStepIndex, KratosCoreFastSuite)
{
    NonHistoricalDataPolicy policy;

    KRATOS_EXPECT_EQ(policy.GetTotalNumberOfSteps(), 1);
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 0);
    KRATOS_EXPECT_TRUE(policy.GetStepCategory() == StepCategory::AnyStep);

    // The single step slot is only reachable through AnyStep
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::AnyStep), 0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        policy.GetStepIndex(StepCategory::TimeStep),
        "NonHistoricalDataPolicy does not support step categories.");

    // Cloning never advances anything
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::TimeStep), 0);
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::AnyStep), 0);
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 0);

    // Clone produces an equal policy
    auto p_clone = policy.Clone();
    KRATOS_EXPECT_TRUE(p_clone->IsSame(policy));
}

KRATOS_TEST_CASE_IN_SUITE(HistoricalDataPolicyStepCycling, KratosCoreFastSuite)
{
    HistoricalDataPolicy policy(StepCategory::TimeStep, 3);

    KRATOS_EXPECT_EQ(policy.GetTotalNumberOfSteps(), 3);
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 0);
    KRATOS_EXPECT_TRUE(policy.GetStepCategory() == StepCategory::TimeStep);

    // Cloning a non-matching category does not advance the ring
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::IterationStep), 0);
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::SubStep), 0);
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 0);

    // Wrap-around indexing with latest index 0: previous steps wrap backwards
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 0), 0);
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 1), 2);
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 2), 1);

    // Cloning the matching category cycles 0 -> 1 -> 2 -> 0
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::TimeStep), 1);
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 1);
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 0), 1);
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 1), 0);
    KRATOS_EXPECT_EQ(policy.GetStepIndex(StepCategory::TimeStep, 2), 2);

    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::TimeStep), 2);
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::TimeStep), 0); // wrap around
    KRATOS_EXPECT_EQ(policy.GetLatestStepIndex(), 0);

    // Clone re-creates the policy configuration but resets the running step index
    KRATOS_EXPECT_EQ(policy.CloneStepIndex(StepCategory::TimeStep), 1);
    auto p_clone = policy.Clone();
    KRATOS_EXPECT_TRUE(p_clone->IsSame(policy));
    KRATOS_EXPECT_EQ(p_clone->GetLatestStepIndex(), 0);
}

} // namespace Testing
} // namespace Kratos
