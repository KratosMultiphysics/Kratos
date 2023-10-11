// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <memory>
#include "strategy_wrapper.hpp"

#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class SolvingStrategyWrapper : public StrategyWrapper
{
public:
    explicit SolvingStrategyWrapper(std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> strategy) :
        mpStrategy(std::move(strategy))
    {
    }

    ~SolvingStrategyWrapper() override = default;

    TimeStepEndState::ConvergenceState GetConvergenceState() override
    {
        return mpStrategy->IsConverged() ? TimeStepEndState::ConvergenceState::converged :
                                           TimeStepEndState::ConvergenceState::non_converged;
    }

    size_t GetNumberOfIterations() const override
    {
        return mpStrategy->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER];
    }

    double GetEndTime() const override
    {
        return mpStrategy->GetModelPart().GetProcessInfo()[TIME];
    }

    void Initialize() override
    {
        mpStrategy->Initialize();
    }

    void InitializeSolutionStep() override
    {
        mpStrategy->InitializeSolutionStep();
    }

    void Predict() override
    {
        mpStrategy->Predict();
    }

    bool SolveSolutionStep() override
    {
        return mpStrategy->SolveSolutionStep();
    }

    void FinalizeSolutionStep() override
    {
        return mpStrategy->FinalizeSolutionStep();
    }

private:
    std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> mpStrategy;
};

}
