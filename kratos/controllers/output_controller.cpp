//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"

// Include base h
#include "output_controller.h"

namespace Kratos {

Parameters OutputController::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name"     : "PLEASE_SPECIFY_THE_MODEL_PART_NAME",
        "output_control_type" : "step",
        "output_interval"     : 1.0
    })" );
}

OutputController::OutputController(
    const Model& rModel,
    Parameters Settings)
    : mpModel(&rModel)
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());

    mModelPartName = Settings["model_part_name"].GetString();

    const auto& output_control_type = Settings["output_control_type"].GetString();
    if (output_control_type == "step") {
        mInterval = Settings["output_interval"].GetInt();
        mpVariable = &STEP;
        mNextPossibleEvaluate = 0;
    } else if (output_control_type == "time") {
        mInterval = Settings["output_interval"].GetDouble();
        mpVariable = &TIME;
        mNextPossibleEvaluate = 0.0;
    } else {
        KRATOS_ERROR << "Unsupported output_control_type = \""
                     << output_control_type << "\". Followings are supported:"
                     << "\n\tstep"
                     << "\n\ttime";
    }

    Update();

    KRATOS_CATCH("");
}

Controller::Pointer OutputController::Create(
    Model& rModel,
    Parameters Settings) const
{
    return Kratos::make_shared<OutputController>(rModel, Settings);
}

int OutputController::Check() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(mModelPartName)) << mModelPartName << " not found in the model.";

    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    std::visit([&](const auto& pVariable){
        KRATOS_ERROR_IF_NOT(r_process_info.Has(*pVariable))
            << pVariable->Name() << " not found in process info of "
            << this->mModelPartName << ".\n";
    }, mpVariable);

    return 0;

    KRATOS_CATCH("");
}

bool OutputController::Evaluate() const
{
    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    return std::visit([&](const auto& pVariable){
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        const auto current_value = r_process_info[*pVariable];
        return std::get<data_type>(mNextPossibleEvaluate) <= current_value;
    }, mpVariable);
}

void OutputController::Update()
{
    KRATOS_TRY

    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    std::visit([&](const auto& pVariable) {
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        const auto interval = std::get<data_type>(mInterval);

        if (interval > 0) {
            const auto current_value = r_process_info[*pVariable];
            auto& next_possible_evaluate = std::get<data_type>(mNextPossibleEvaluate);
            next_possible_evaluate = (static_cast<int>(current_value / interval) + 1) * interval;
        }
    }, mpVariable);

    KRATOS_CATCH("");
}

std::variant<int, double> OutputController::GetCurrentControlValue() const
{
    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    return std::visit([&](const auto& pVariable) -> std::variant<int, double> {
        return r_process_info[*pVariable];
    }, mpVariable);
}

std::variant<int, double> OutputController::GetInterval() const
{
    return mInterval;
}

std::variant<int, double> OutputController::GetNextPossibleEvaluateControlValue() const
{
    return mNextPossibleEvaluate;
}

std::string OutputController::Info() const
{
    return "OutputController";
}

void OutputController::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void OutputController::PrintData(std::ostream& rOStream) const
{
    std::visit([&](const auto& pVariable){
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        rOStream << "Interval     : " << std::get<data_type>(mInterval) << std::endl
                 << "Variable     : " << pVariable->Name() << std::endl
                 << "Next evaluate: " << std::get<data_type>(mNextPossibleEvaluate) << std::endl;
    }, mpVariable);
}

} // namespace Kratos.
