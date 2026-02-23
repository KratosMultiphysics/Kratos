//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jonathan Nuttall
//
#include "set_multiple_moving_loads.h"
#include "includes/variables.h"

namespace Kratos
{
SetMultipleMovingLoadsProcess::SetMultipleMovingLoadsProcess(ModelPart& rModelPart, const Parameters& rProcessSettings)
    : mrModelPart(rModelPart), mParameters(rProcessSettings)
{
    Parameters default_parameters(R"(
        {
            "help"                    : "This process applies multiple coupled moving load conditions belonging to a model part. The loads move over line elements.",
            "model_part_name"         : "please_specify_model_part_name",
            "compute_model_part_name" : "please_specify_compute_body_part_name",
            "variable_name"           : "POINT_LOAD",
            "load"                    : [0.0, 1.0, 0.0],
            "direction"               : [1,1,1],
            "velocity"                : 1,
            "origin"                  : [0.0, 0.0, 0.0],
            "configuration"           : [0.0]
        }  )");

    // Set default velocity as a string, if the input velocity is a string
    if (mParameters.Has("velocity") && mParameters["velocity"].IsString()) {
        default_parameters["velocity"].SetString("1");
    }

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // check if load parameter has size 3
    KRATOS_ERROR_IF(mParameters["load"].size() != 3) << "'load' has to have size 3!" << std::endl;

    // check if all elements in load parameter are either string or a number
    bool is_all_string = true;
    bool is_all_number = true;
    for (IndexType i = 0; i < mParameters["load"].size(); i++) {
        if (!mParameters["load"][i].IsString()) {
            is_all_string = false;
        }
        if (!mParameters["load"][i].IsNumber()) {
            is_all_number = false;
        }
    }

    KRATOS_ERROR_IF(!is_all_string && !is_all_number)
        << "'load' has to be a vector of numbers, or an array with strings" << std::endl;

    int count = 0;
    for (double offset : mParameters["configuration"].GetVector()) {
        auto parameters_moving_load = mParameters.Clone();

        count++;
        const std::string& newModelPartName = mrModelPart.Name() + "_cloned_" + std::to_string(count);
        auto& new_cloned_model_part = CloneMovingConditionInComputeModelPart(newModelPartName);

        parameters_moving_load.RemoveValue("configuration");
        parameters_moving_load.RemoveValue("compute_model_part_name");
        parameters_moving_load.AddDouble("offset", offset);

        auto moving_point_process =
            Kratos::make_unique<SetMovingLoadProcess>(new_cloned_model_part, parameters_moving_load);

        mMovingPointLoadsProcesses.push_back(std::move(moving_point_process));
    }

    RemoveClonedConditions();
}

ModelPart& SetMultipleMovingLoadsProcess::CloneMovingConditionInComputeModelPart(const std::string& NewBodyPartName)
{
    auto& compute_model_part =
        mrModelPart.GetRootModelPart().GetSubModelPart(mParameters["compute_model_part_name"].GetString());
    auto& new_model_part = compute_model_part.CreateSubModelPart(NewBodyPartName);
    new_model_part.AddNodes(mrModelPart.NodesBegin(), mrModelPart.NodesEnd());

    int index = GetMaxConditionsIndex();
    for (auto& moving_load_condition : mrModelPart.Conditions()) {
        index++;
        const Condition::Pointer CloneCondition =
            moving_load_condition.Clone(index, moving_load_condition.GetGeometry());
        new_model_part.AddCondition(CloneCondition);
    }
    return new_model_part;
}

int SetMultipleMovingLoadsProcess::GetMaxConditionsIndex() const
{
    int max_index = 0;
    for (const auto& condition : mrModelPart.GetRootModelPart().Conditions()) {
        max_index = std::max(max_index, static_cast<int>(condition.Id()));
    }
    return max_index;
}

void SetMultipleMovingLoadsProcess::RemoveClonedConditions()
{
    auto& compute_model_part =
        mrModelPart.GetRootModelPart().GetSubModelPart(mParameters["compute_model_part_name"].GetString());

    for (const auto& moving_load_condition : mrModelPart.Conditions()) {
        if (compute_model_part.HasCondition(moving_load_condition))
            compute_model_part.pGetCondition(moving_load_condition)->Set(TO_ERASE, true);
    }
    // Call method
    compute_model_part.RemoveConditions(TO_ERASE);
}

void SetMultipleMovingLoadsProcess::ExecuteInitialize()
{
    for (const auto& rMovingPointLoad : mMovingPointLoadsProcesses) {
        rMovingPointLoad->ExecuteInitialize();
    }
}

void SetMultipleMovingLoadsProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& rMovingPointLoad : mMovingPointLoadsProcesses) {
        rMovingPointLoad->ExecuteInitializeSolutionStep();
    }
}

void SetMultipleMovingLoadsProcess::ExecuteFinalizeSolutionStep()
{
    for (const auto& rMovingPointLoad : mMovingPointLoadsProcesses) {
        rMovingPointLoad->ExecuteFinalizeSolutionStep();
    }
}

std::string SetMultipleMovingLoadsProcess::Info() const { return "SetMultipleMovingLoadsProcess"; }
} // namespace Kratos