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
#include "solving_strategies/strategies/solving_strategy.h"

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

class NeverConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
   bool IsConverged() override {return false;}
};

class AlwaysConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
    bool IsConverged() override {return true;}
};

}

class StubFixedSizeStepIncrementor : public TimeIncrementor
{
public:
    double GetIncrement() override {return 1.23;}
};

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsStartStateandPerformedStatesWhenNoStartStateIsGiven, KratosGeoMechanicsWip)
{
    TimeLoopExecutor<AlwaysConvergingSolverStrategy> executor;
    executor.SetTimeIncrementor(std::make_unique<StubFixedSizesStepIncrementor>());
    const auto step_states = executor.Run();

    KRATOS_EXPECT_EQ(1, step_states.size());
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsOnlyPerformedStatesWhenStartStateIsGiven, KratosGeoMechanicsWip)
{
    TimeLoopExecutor<AlwaysConvergingSolverStrategy> executor;
    executor.SetTimeIncrementor(std::make_unique<StubFixedSizesStepIncrementor>());
    const auto step_states = executor.Run();

    KRATOS_EXPECT_EQ(1, step_states.size());
}

}