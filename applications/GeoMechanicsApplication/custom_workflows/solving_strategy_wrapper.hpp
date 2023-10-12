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

    void SetEndTime(double EndTime) override
    {

    }

    double GetTimeIncrement() const override
    {
        return 0;
    }

    void SetTimeIncrement(double TimeIncrement) override
    {

    }

    size_t GetStepNumber() const override
    {
        return 0;
    }

    void IncrementStepNumber() override
    {

    }

    void CloneTimeStep() override
    {

    }

    void RestorePositionsAndDOFVectorToStartOfStep() override
    {

    }

    void SaveTotalDisplacementFieldAtStartOfStage() override
    {

    }

    void AccumulateTotalDisplacementField() override
    {

    }

    void OutputProcess() override
    {

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
