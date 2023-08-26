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
#include "temporal_controller.h"

namespace Kratos {

Parameters TemporalController::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name"     : "PLEASE_SPECIFY_THE_MODEL_PART_NAME",
        "output_control_type" : "step",
        "output_interval"     : 1.0
    })" );
}

TemporalController::TemporalController(
    const Model& rModel,
    Parameters Settings)
    : mpModel(&rModel),
      mNextOutput{}
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());

    mModelPartName = Settings["model_part_name"].GetString();

    const auto& output_control_type = Settings["output_control_type"].GetString();
    if (output_control_type == "step") {
        mInterval = Settings["output_interval"].GetInt();
        mpVariable = &STEP;
    } else if (output_control_type == "time") {
        mInterval = Settings["output_interval"].GetDouble();
        mpVariable = &TIME;
    } else {
        KRATOS_ERROR << "Unsupported output_control_type = \""
                     << output_control_type << "\". Followings are supported:"
                     << "\n\tstep"
                     << "\n\ttime";
    }

    ScheduleNextOutput();

    KRATOS_CATCH("");
}

void TemporalController::ScheduleNextOutput()
{
    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    std::visit([&](const auto& pVariable) {
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        const auto interval = std::get<data_type>(mInterval);

        if (interval > 0) {
            const auto current_value = r_process_info[*pVariable];
            while (std::get<data_type>(mNextOutput) <= current_value) {
                std::get<data_type>(mNextOutput) += interval;
            }
        }
    }, mpVariable);
}

Controller::Pointer TemporalController::Create(
    Model& rModel,
    Parameters Settings) const
{
    return Kratos::make_shared<TemporalController>(rModel, Settings);
}

int TemporalController::Check() const
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

bool TemporalController::Evaluate()
{
    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    return std::visit([&](const auto& pVariable){
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        const auto current_value = r_process_info[*pVariable];
        if (std::get<data_type>(mNextOutput) <= current_value) {
            ScheduleNextOutput();
            return true;
        } else {
            return false;
        }
    }, mpVariable);
}

std::variant<int, double> TemporalController::GetCurrentControlValue() const
{
    const auto& r_process_info = mpModel->GetModelPart(mModelPartName).GetProcessInfo();
    return std::visit([&](const auto& pVariable) -> std::variant<int, double> {
        return r_process_info[*pVariable];
    }, mpVariable);
}

std::string TemporalController::Info() const
{
    return "TemporalController";
}

void TemporalController::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void TemporalController::PrintData(std::ostream& rOStream) const
{
    std::visit([&](const auto& pVariable){
        using data_type = typename std::decay_t<decltype(*pVariable)>::Type;
        rOStream << "Interval   : " << std::get<data_type>(mInterval) << std::endl
                 << "Variable   : " << pVariable->Name() << std::endl
                 << "Next output: " << std::get<data_type>(mNextOutput) << std::endl;
    }, mpVariable);
}

} // namespace Kratos.
