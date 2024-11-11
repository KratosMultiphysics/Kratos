// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_workflows/strategy_wrapper.hpp"
#include "custom_workflows/time_step_executor.h"
#include "geo_mechanics_fast_suite.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include <gmock/gmock.h>

using namespace Kratos;

namespace
{

class ProcessSpy : public Process
{
public:
    MOCK_METHOD(void, ExecuteInitializeSolutionStep, (), (override));
    MOCK_METHOD(void, ExecuteFinalizeSolutionStep, (), (override));
};

class DummyStrategyWrapper : public StrategyWrapper
{
public:
    explicit DummyStrategyWrapper(TimeStepEndState::ConvergenceState ConvergenceState)
        : mConvergenceState(ConvergenceState)
    {
    }

    MOCK_METHOD(std::size_t, GetNumberOfIterations, (), (const, override));

    [[nodiscard]] double GetEndTime() const override { return 10.; };

    MOCK_METHOD(void, SetEndTime, (double EndTime), (override));

    [[nodiscard]] double GetTimeIncrement() const override { return 0.0; }

    MOCK_METHOD(void, SetTimeIncrement, (double TimeIncrement), (override));

    [[nodiscard]] std::size_t GetStepNumber() const override { return 0; }

    MOCK_METHOD(void, IncrementStepNumber, (), (override));

    MOCK_METHOD(void, CloneTimeStep, (), (override));

    MOCK_METHOD(void, RestorePositionsAndDOFVectorToStartOfStep, (), (override));

    MOCK_METHOD(void, SaveTotalDisplacementFieldAtStartOfTimeLoop, (), (override));

    MOCK_METHOD(void, AccumulateTotalDisplacementField, (), (override));

    MOCK_METHOD(void, ComputeIncrementalDisplacementField, (),  (override));

    MOCK_METHOD(void, OutputProcess, (), (override));

    void Initialize() override { ++mSolverStrategyInitializeCalls; }

    MOCK_METHOD(void, InitializeOutput, (), (override));

    void InitializeSolutionStep() override { ++mSolverStrategyInitializeSolutionStepCalls; }

    void Predict() override { ++mSolverStrategyPredictCalls; }

    TimeStepEndState::ConvergenceState SolveSolutionStep() override
    {
        ++mSolverStrategySolveSolutionsStepCalls;
        return mConvergenceState;
    }

    void FinalizeSolutionStep() override { ++mSolverStrategyFinalizeSolutionStepCalls; }

    [[nodiscard]] unsigned int NumberOfSolverStrategyInitializeCalls() const
    {
        return mSolverStrategyInitializeCalls;
    }

    [[nodiscard]] unsigned int NumberOfSolverStrategyInitializeSolutionStepCalls() const
    {
        return mSolverStrategyInitializeSolutionStepCalls;
    }

    [[nodiscard]] unsigned int NumberOfSolverStrategyPredictCalls() const
    {
        return mSolverStrategyPredictCalls;
    }

    [[nodiscard]] unsigned int NumberOfSolverStrategySolveSolutionStepCalls() const
    {
        return mSolverStrategySolveSolutionsStepCalls;
    }

    [[nodiscard]] unsigned int NumberOfSolverStrategyFinalizeSolutionStepCalls() const
    {
        return mSolverStrategyFinalizeSolutionStepCalls;
    }

    MOCK_METHOD(void, FinalizeOutput, (), (override));

private:
    TimeStepEndState::ConvergenceState mConvergenceState;
    unsigned int                       mSolverStrategyInitializeCalls             = 0;
    unsigned int                       mSolverStrategyInitializeSolutionStepCalls = 0;
    unsigned int                       mSolverStrategyPredictCalls                = 0;
    unsigned int                       mSolverStrategySolveSolutionsStepCalls     = 0;
    unsigned int                       mSolverStrategyFinalizeSolutionStepCalls   = 0;
};

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RunReturnsNonConvergedWhenStrategyDoesNotConverge, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             nonconverging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::non_converged);
    executor.SetSolverStrategy(nonconverging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_TRUE(executor.Run(time).NonConverged())
}

KRATOS_TEST_CASE_IN_SUITE(RunReturnsConvergedWhenStrategyConverged, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_TRUE(executor.Run(time).Converged())
}

KRATOS_TEST_CASE_IN_SUITE(ProcessMemberFunctionsAllCalledOnce, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    auto spy = std::make_shared<ProcessSpy>();
    EXPECT_CALL(*spy, ExecuteInitializeSolutionStep()).Times(1);
    EXPECT_CALL(*spy, ExecuteFinalizeSolutionStep()).Times(1);

    std::vector<std::shared_ptr<Process>> processes{spy};
    std::vector<std::weak_ptr<Process>>   process_observables{spy};
    executor.SetProcessObservables(process_observables);
    const auto time = 0.0;

    executor.Run(time);
}

KRATOS_TEST_CASE_IN_SUITE(SolverStrategyMemberFunctionsAllExceptInitializeAndFinalizeCalledOnce,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;

    executor.Run(time);

    KRATOS_EXPECT_EQ(0, converging_strategy->NumberOfSolverStrategyInitializeCalls());
    KRATOS_EXPECT_EQ(1, converging_strategy->NumberOfSolverStrategyInitializeSolutionStepCalls());
    KRATOS_EXPECT_EQ(1, converging_strategy->NumberOfSolverStrategyPredictCalls());
    KRATOS_EXPECT_EQ(1, converging_strategy->NumberOfSolverStrategySolveSolutionStepCalls());
    KRATOS_EXPECT_EQ(0, converging_strategy->NumberOfSolverStrategyFinalizeSolutionStepCalls());
}

KRATOS_TEST_CASE_IN_SUITE(ConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(NonConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             non_converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::non_converged);
    executor.SetSolverStrategy(non_converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(TimeStepExecutionReturnsNumberOfIterations, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeStepExecutor executor;
    auto             converging_strategy =
        std::make_shared<DummyStrategyWrapper>(TimeStepEndState::ConvergenceState::converged);
    EXPECT_CALL(*converging_strategy, GetNumberOfIterations()).WillOnce(testing::Return(4));

    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_EQ(4, executor.Run(time).num_of_iterations);
}

} // namespace Kratos::Testing
