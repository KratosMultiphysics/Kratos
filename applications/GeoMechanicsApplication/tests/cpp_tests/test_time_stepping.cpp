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
#include "testing/testing.h"

using namespace Kratos;

namespace
{

class ProcessSpy : public Process
{
public:
    void ExecuteInitializeSolutionStep() override
    {
        ++mSolutionStepInitializedCalls;
    }
    void ExecuteFinalizeSolutionStep() override
    {
        ++mSolutionStepFinalizedCalls;
    }

    [[nodiscard]] unsigned int NumberOfExecuteInitializeSolutionStepCalls() const
    {
        return mSolutionStepInitializedCalls;
    }
    [[nodiscard]] unsigned int NumberOfExecuteFinalizeSolutionStepCalls() const
    {
        return mSolutionStepFinalizedCalls;
    }

private:
    unsigned int mSolutionStepInitializedCalls = 0;
    unsigned int mSolutionStepFinalizedCalls = 0;
};

class DummyStrategyWrapper : public StrategyWrapper
{
public:
    explicit DummyStrategyWrapper(TimeStepEndState::ConvergenceState ConvergenceState)
        : mConvergenceState(ConvergenceState)
    {
    }
    [[nodiscard]] std::size_t GetNumberOfIterations() const override
    {
        return 4;
    };
    [[nodiscard]] double GetEndTime() const override { return 10.; };
    void SetEndTime(double EndTime) override
    {
        // intentionally empty
    }
    [[nodiscard]] double GetTimeIncrement() const override { return 0.0; }
    void SetTimeIncrement(double TimeIncrement) override
    {
        // intentionally empty
    }
    [[nodiscard]] std::size_t GetStepNumber() const override { return 0; }
    void IncrementStepNumber() override
    {
        // intentionally empty
    }
    void CloneTimeStep() override{
        // intentionally empty
    };
    void RestorePositionsAndDOFVectorToStartOfStep() override{
        // intentionally empty
    };
    void SaveTotalDisplacementFieldAtStartOfTimeLoop() override{
        // intentionally empty
    };
    void AccumulateTotalDisplacementField() override{
        // intentionally empty
    };
    void OutputProcess() override{
        // intentionally empty
    };

    void Initialize() override { ++mSolverStrategyInitializeCalls; }

    void InitializeOutput() override
    {
        // Intentionally empty
    }

    void InitializeSolutionStep() override
    {
        ++mSolverStrategyInitializeSolutionStepCalls;
    }
    void Predict() override { ++mSolverStrategyPredictCalls; }
    TimeStepEndState::ConvergenceState SolveSolutionStep() override
    {
        ++mSolverStrategySolveSolutionsStepCalls;
        return mConvergenceState;
    }
    void FinalizeSolutionStep() override
    {
        ++mSolverStrategyFinalizeSolutionStepCalls;
    }
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
    void FinalizeOutput() override
    {
        // intentionally empty
    }

private:
    TimeStepEndState::ConvergenceState mConvergenceState;
    unsigned int mSolverStrategyInitializeCalls = 0;
    unsigned int mSolverStrategyInitializeSolutionStepCalls = 0;
    unsigned int mSolverStrategyPredictCalls = 0;
    unsigned int mSolverStrategySolveSolutionsStepCalls = 0;
    unsigned int mSolverStrategyFinalizeSolutionStepCalls = 0;
};

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RunReturnsNonConvergedWhenStrategyDoesNotConverge, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto nonconverging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::non_converged);
    executor.SetSolverStrategy(nonconverging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_TRUE(executor.Run(time).NonConverged())
}

KRATOS_TEST_CASE_IN_SUITE(RunReturnsConvergedWhenStrategyConverged, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_TRUE(executor.Run(time).Converged())
}

KRATOS_TEST_CASE_IN_SUITE(ProcessMemberFunctionsAllCalledOnce, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    auto spy = std::make_shared<ProcessSpy>();

    std::vector<std::shared_ptr<Process>> processes{spy};
    std::vector<std::weak_ptr<Process>> process_observables{spy};
    executor.SetProcessObservables(process_observables);
    const auto time = 0.0;

    executor.Run(time);

    KRATOS_EXPECT_EQ(1, spy->NumberOfExecuteInitializeSolutionStepCalls());
    KRATOS_EXPECT_EQ(1, spy->NumberOfExecuteFinalizeSolutionStepCalls());
}

KRATOS_TEST_CASE_IN_SUITE(SolverStrategyMemberFunctionsAllExceptInitializeAndFinalizeCalledOnce,
                          KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;

    executor.Run(time);

    KRATOS_EXPECT_EQ(0, converging_strategy->NumberOfSolverStrategyInitializeCalls());
    KRATOS_EXPECT_EQ(
        1, converging_strategy->NumberOfSolverStrategyInitializeSolutionStepCalls());
    KRATOS_EXPECT_EQ(1, converging_strategy->NumberOfSolverStrategyPredictCalls());
    KRATOS_EXPECT_EQ(1, converging_strategy->NumberOfSolverStrategySolveSolutionStepCalls());
    KRATOS_EXPECT_EQ(
        0, converging_strategy->NumberOfSolverStrategyFinalizeSolutionStepCalls());
}

KRATOS_TEST_CASE_IN_SUITE(ConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(NonConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto non_converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::non_converged);
    executor.SetSolverStrategy(non_converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(TimeStepExecutionReturnsNumberOfIterations, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyStrategyWrapper>(
        TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_EQ(4, executor.Run(time).num_of_iterations);
}

} // namespace Kratos::Testing
