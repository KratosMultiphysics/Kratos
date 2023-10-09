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

#include "testing/testing.h"
#include "custom_workflows/time_loop_executor.hpp"
#include "custom_workflows/time_step_end_state.hpp"
#include "custom_workflows/prescribed_time_incrementor.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include <numeric>

using namespace Kratos;

namespace
{

class DummySparseType
{
public:
    struct DataType{};
    struct MatrixType{};
    struct VectorType{};
    struct MatrixPointerType{};
    struct VectorPointerType{};
};

class DummyDenseType
{
public:
    struct MatrixType{};
    struct VectorType{};
};

class AlwaysConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
    bool IsConverged() override {return true;}
};

class DummySolverStrategy : public StrategyWrapper
{
public:
    DummySolverStrategy(TimeStepEndState::ConvergenceState AConvergenceState) : mConvergenceState{AConvergenceState} {}
    [[nodiscard]] TimeStepEndState::ConvergenceState GetConvergenceState()         override
    {
        return mConvergenceState;
    }
    [[nodiscard]] std::size_t                        GetNumberOfIterations() const override
    {
        return 1;
    }
    [[nodiscard]] double                             GetEndTime() const            override
    {
        return 0.0;
    }
    void Initialize()             override {}
    void InitializeSolutionStep() override {}
    void Predict()                override {}
    bool SolveSolutionStep()      override { return false; }
    void FinalizeSolutionStep()   override {}

private:
    TimeStepEndState::ConvergenceState mConvergenceState{TimeStepEndState::ConvergenceState::converged};
};
}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsPerformedStatesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_EQ(increments.size(), step_states.size());
    KRATOS_EXPECT_TRUE(std::all_of(step_states.begin(), step_states.end(), [](const auto& step_state) { return step_state.Converged(); }))
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsOneNonConvergedPerformedStatesAfterRunningAnNeverConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::non_converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_EQ(1, step_states.size());
    KRATOS_EXPECT_TRUE(step_states[0].NonConverged());
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsEndTimesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_DOUBLE_EQ(std::accumulate(increments.begin(), increments.end(), 0.), step_states.back().time);
}

}