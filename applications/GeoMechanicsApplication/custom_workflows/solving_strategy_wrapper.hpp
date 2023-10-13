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
#include "geo_mechanics_application_variables.h"
#include "write_output.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class SolvingStrategyWrapper : public StrategyWrapper
{
public:
    SolvingStrategyWrapper(std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> strategy,
                           bool ResetDisplacements = false,
                           const std::filesystem::path& rWorkingDirectory = "",
                           const Parameters& rProjectParameters = {}) :
        mpStrategy(std::move(strategy)),
        mrModelPart(mpStrategy->GetModelPart()),
        mResetDisplacements{ResetDisplacements},
        mProjectParameters{rProjectParameters},
        mWorkingDirectory{rWorkingDirectory}
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
        return mrModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
    }

    double GetEndTime() const override
    {
        return mrModelPart.GetProcessInfo()[TIME];
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
        mrModelPart.GetProcessInfo()[TIME] = EndTime;
    }

    double GetTimeIncrement() const override
    {
        return mrModelPart.GetProcessInfo()[DELTA_TIME];
    }

    void SetTimeIncrement(double TimeIncrement) override
    {
        mrModelPart.GetProcessInfo()[DELTA_TIME] = TimeIncrement;
    }

    size_t GetStepNumber() const override
    {
        return static_cast<std::size_t>(mrModelPart.GetProcessInfo()[STEP]);
    }

    void IncrementStepNumber() override
    {
        mrModelPart.GetProcessInfo()[STEP] += 1;
    }

    void CloneTimeStep() override
    {
        mrModelPart.CloneTimeStep();
    }

    void RestorePositionsAndDOFVectorToStartOfStep() override
    {
        VariableUtils().UpdateCurrentPosition(mrModelPart.Nodes(), DISPLACEMENT, 1);
        for (auto& node : mrModelPart.Nodes())
        {
            node.GetSolutionStepValue(DISPLACEMENT, 0) = node.GetSolutionStepValue(DISPLACEMENT, 1);
            node.GetSolutionStepValue(WATER_PRESSURE, 0) = node.GetSolutionStepValue(WATER_PRESSURE, 1);
            if (node.Has(ROTATION)) node.GetSolutionStepValue(ROTATION, 0) = node.GetSolutionStepValue(ROTATION, 1);
        }
    }

    void SaveTotalDisplacementFieldAtStartOfStage() override
    {
        if (mResetDisplacements)
        {
            for (const auto& node : mrModelPart.Nodes())
            {
                mOldTotalDisplacements.emplace_back(node.GetSolutionStepValue(TOTAL_DISPLACEMENT));
            }
        }
    }

    void AccumulateTotalDisplacementField() override
    {
        if (mResetDisplacements)
        {
            KRATOS_ERROR_IF_NOT(mrModelPart.Nodes().size() == mOldTotalDisplacements.size()) << "TEST";
            std::size_t count = 0;
            for (auto& node : mrModelPart.Nodes())
            {
                node.GetSolutionStepValue(TOTAL_DISPLACEMENT) = mOldTotalDisplacements[count] + node.GetSolutionStepValue(DISPLACEMENT);
                count++;
            }
        }

    }

    void OutputProcess() override
    {
        GeoOutputWriter::WriteGiDOutput(mrModelPart, mProjectParameters, mWorkingDirectory.generic_string());
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
    ModelPart& mrModelPart;
    bool mResetDisplacements;
    std::vector<array_1d<double, 3>> mOldTotalDisplacements;
    std::filesystem::path mWorkingDirectory;
    Parameters mProjectParameters;
};

}
