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
#include "solving_strategies/strategies/solving_strategy.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

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
    DummyStrategyWrapper()
    {
        ON_CALL(*this, GetNumberOfIterations()).WillByDefault(testing::Return(4));
        ON_CALL(*this, GetEndTime()).WillByDefault(testing::Return(10.0));
        ON_CALL(*this, GetTimeIncrement()).WillByDefault(testing::Return(0.0));
        ON_CALL(*this, GetStepNumber()).WillByDefault(testing::Return(0));
    }

    MOCK_METHOD(std::size_t, GetNumberOfIterations, (), (const, override));
    MOCK_METHOD(double, GetEndTime, (), (const, override));
    MOCK_METHOD(void, SetEndTime, (double EndTime), (override));
    MOCK_METHOD(double, GetTimeIncrement, (), (const, override));
    MOCK_METHOD(void, SetTimeIncrement, (double TimeIncrement), (override));
    MOCK_METHOD(std::size_t, GetStepNumber, (), (const, override));
    MOCK_METHOD(void, IncrementStepNumber, (), (override));
    MOCK_METHOD(void, CloneTimeStep, (), (override));
    MOCK_METHOD(void, RestorePositionsAndDOFVectorToStartOfStep, (), (override));
    MOCK_METHOD(void, AccumulateTotalDisplacementField, (), (override));
    MOCK_METHOD(void, ComputeIncrementalDisplacementField, (), (override));
    MOCK_METHOD(void, OutputProcess, (), (override));
    MOCK_METHOD(void, Initialize, (), (override));
    MOCK_METHOD(void, InitializeOutput, (), (override));
    MOCK_METHOD(void, InitializeSolutionStep, (), (override));
    MOCK_METHOD(void, Predict, (), (override));
    MOCK_METHOD(TimeStepEndState::ConvergenceState, SolveSolutionStep, (), (override));
    MOCK_METHOD(void, FinalizeSolutionStep, (), (override));
    MOCK_METHOD(void, FinalizeOutput, (), (override));
};

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, RunReturnsNonConvergedWhenStrategyDoesNotConverge)
{
    TimeStepExecutor executor;
    auto             nonconverging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*nonconverging_strategy, SolveSolutionStep())
        .WillOnce(testing::Return(TimeStepEndState::ConvergenceState::non_converged));

    executor.SetSolverStrategy(nonconverging_strategy);
    constexpr auto time = 0.0;
    EXPECT_TRUE(executor.Run(time).NonConverged());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, RunReturnsConvergedWhenStrategyConverged)
{
    TimeStepExecutor executor;
    auto             converging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::converged));

    executor.SetSolverStrategy(converging_strategy);
    constexpr auto time = 0.0;
    EXPECT_TRUE(executor.Run(time).Converged());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ProcessMemberFunctionsAllCalledOnce)
{
    TimeStepExecutor executor;
    auto             converging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::converged));

    executor.SetSolverStrategy(converging_strategy);
    auto spy = std::make_shared<ProcessSpy>();
    EXPECT_CALL(*spy, ExecuteInitializeSolutionStep()).Times(1);
    EXPECT_CALL(*spy, ExecuteFinalizeSolutionStep()).Times(1);

    std::vector<std::shared_ptr<Process>> processes{spy};
    std::vector<std::weak_ptr<Process>>   process_observables{spy};
    executor.SetProcessObservables(process_observables);
    constexpr auto time = 0.0;

    executor.Run(time);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, SolverStrategyMemberFunctionsAllExceptInitializeAndFinalizeCalledOnce)
{
    TimeStepExecutor executor;
    auto             converging_strategy = std::make_shared<DummyStrategyWrapper>();

    executor.SetSolverStrategy(converging_strategy);
    constexpr auto time = 0.0;

    EXPECT_CALL(*converging_strategy, Initialize()).Times(0);
    EXPECT_CALL(*converging_strategy, InitializeSolutionStep()).Times(1);
    EXPECT_CALL(*converging_strategy, Predict()).Times(1);
    EXPECT_CALL(*converging_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::converged));
    EXPECT_CALL(*converging_strategy, FinalizeSolutionStep()).Times(0);

    executor.Run(time);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ConvergingTimeStepExecutionReturnsGivenTime)
{
    TimeStepExecutor executor;
    auto             converging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::converged));

    executor.SetSolverStrategy(converging_strategy);
    constexpr auto time = 2.0;
    EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NonConvergingTimeStepExecutionReturnsGivenTime)
{
    TimeStepExecutor executor;
    auto             non_converging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*non_converging_strategy, SolveSolutionStep())
        .WillOnce(testing::Return(TimeStepEndState::ConvergenceState::non_converged));

    executor.SetSolverStrategy(non_converging_strategy);
    constexpr auto time = 2.0;
    EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TimeStepExecutionReturnsNumberOfIterations)
{
    TimeStepExecutor executor;
    auto             converging_strategy = std::make_shared<DummyStrategyWrapper>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::converged));

    EXPECT_CALL(*converging_strategy, GetNumberOfIterations()).WillOnce(testing::Return(4));

    executor.SetSolverStrategy(converging_strategy);
    constexpr auto time = 2.0;
    EXPECT_EQ(4, executor.Run(time).num_of_iterations);
}

} // namespace Kratos::Testing
