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

#include "strategy_wrapper.hpp"
#include <memory>

#include "geo_mechanics_application_variables.h"
#include "geo_output_writer.h"
#include "includes/variables.h"
#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class SolvingStrategyWrapper : public StrategyWrapper
{
public:
    explicit SolvingStrategyWrapper(std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> strategy,
                                    bool                         ResetDisplacements = false,
                                    const std::filesystem::path& rWorkingDirectory  = "",
                                    const Parameters&            rProjectParameters = {})
        : mpStrategy(std::move(strategy)),
          mrModelPart(mpStrategy->GetModelPart()),
          mResetDisplacements{ResetDisplacements},
          mProjectParameters{rProjectParameters},
          mWorkingDirectory{rWorkingDirectory}
    {
    }

    ~SolvingStrategyWrapper() override = default;

    size_t GetNumberOfIterations() const override
    {
        return mrModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
    }

    double GetEndTime() const override { return mrModelPart.GetProcessInfo()[TIME]; }

    void Initialize() override { mpStrategy->Initialize(); }

    void InitializeOutput() override
    {
        const auto gid_output_settings =
            mProjectParameters["output_processes"]["gid_output"][0]["Parameters"];
        mWriter = std::make_unique<GeoOutputWriter>(
            gid_output_settings, mWorkingDirectory.generic_string(), mrModelPart);
    }

    void InitializeSolutionStep() override { mpStrategy->InitializeSolutionStep(); }

    void Predict() override { mpStrategy->Predict(); }

    void SetEndTime(double EndTime) override { mrModelPart.GetProcessInfo()[TIME] = EndTime; }

    double GetTimeIncrement() const override { return mrModelPart.GetProcessInfo()[DELTA_TIME]; }

    void SetTimeIncrement(double TimeIncrement) override
    {
        mrModelPart.GetProcessInfo()[DELTA_TIME] = TimeIncrement;
    }

    size_t GetStepNumber() const override
    {
        return static_cast<std::size_t>(mrModelPart.GetProcessInfo()[STEP]);
    }

    void IncrementStepNumber() override { mrModelPart.GetProcessInfo()[STEP] += 1; }

    void CloneTimeStep() override
    {
        auto& root_model_part = mrModelPart.IsSubModelPart() ? mrModelPart.GetRootModelPart() : mrModelPart;
        root_model_part.CloneTimeStep();
    }

    void RestorePositionsAndDOFVectorToStartOfStep() override
    {
        VariableUtils().UpdateCurrentPosition(mrModelPart.Nodes(), DISPLACEMENT, 1);

        const auto index_of_old_value = Node::IndexType{1};
        const auto index_of_new_value = Node::IndexType{0};
        CopyNodalSolutionStepValues(DISPLACEMENT, index_of_old_value, index_of_new_value);
        CopyNodalSolutionStepValues(WATER_PRESSURE, index_of_old_value, index_of_new_value);
        CopyNodalSolutionStepValues(ROTATION, index_of_old_value, index_of_new_value);
    }

    void SaveTotalDisplacementFieldAtStartOfTimeLoop() override
    {
        if (mResetDisplacements) {
            mOldTotalDisplacements.clear();
            for (const auto& node : mrModelPart.Nodes()) {
                mOldTotalDisplacements.emplace_back(node.GetSolutionStepValue(TOTAL_DISPLACEMENT));
            }
        }
    }

    void AccumulateTotalDisplacementField() override
    {
        if (mResetDisplacements) {
            KRATOS_ERROR_IF_NOT(mrModelPart.Nodes().size() == mOldTotalDisplacements.size())
                << "The number of old displacements (" << mOldTotalDisplacements.size()
                << ") does not match the current number of nodes (" << mrModelPart.Nodes().size() << ").";
            std::size_t count = 0;
            for (auto& node : mrModelPart.Nodes()) {
                node.GetSolutionStepValue(TOTAL_DISPLACEMENT) =
                    mOldTotalDisplacements[count] + node.GetSolutionStepValue(DISPLACEMENT);
                ++count;
            }
        }
    }

    void OutputProcess() override
    {
        if (mWriter) {
            const auto write_hydraulic_head_to_nodes = false;
            const auto gid_output_settings =
                mProjectParameters["output_processes"]["gid_output"][0]["Parameters"];
            mWriter->WriteGiDOutput(mrModelPart, gid_output_settings, write_hydraulic_head_to_nodes);
        }
    }

    TimeStepEndState::ConvergenceState SolveSolutionStep() override
    {
        return mpStrategy->SolveSolutionStep() ? TimeStepEndState::ConvergenceState::converged
                                               : TimeStepEndState::ConvergenceState::non_converged;
        ;
    }

    void FinalizeSolutionStep() override { return mpStrategy->FinalizeSolutionStep(); }

    void FinalizeOutput() override
    {
        if (mWriter) {
            mWriter->FinalizeResults();
        }
    }

private:
    template <typename TVariableType>
    void CopyNodalSolutionStepValues(const TVariableType& rVariable, Node::IndexType SourceIndex, Node::IndexType DestinationIndex)
    {
        if (!mrModelPart.HasNodalSolutionStepVariable(rVariable)) return;

        for (auto& node : mrModelPart.Nodes()) {
            node.GetSolutionStepValue(rVariable, DestinationIndex) =
                node.GetSolutionStepValue(rVariable, SourceIndex);
        }
    }

    std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> mpStrategy;
    ModelPart&                                                  mrModelPart;
    bool                                                        mResetDisplacements;
    std::vector<array_1d<double, 3>>                            mOldTotalDisplacements;
    Parameters                                                  mProjectParameters;
    std::filesystem::path                                       mWorkingDirectory;
    std::unique_ptr<GeoOutputWriter>                            mWriter;
};

} // namespace Kratos
