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
#include <sstream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "rans_compute_previous_step_bossak_second_derivatives_process.h"

namespace Kratos
{
template<class TVariableDataType>
RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::BossakPreviousStepSecondDerivativeCalculator(
    const double BossakAlpha,
    const Variable<TVariableDataType>& rSecondDerivative,
    const Variable<TVariableDataType>& rRelaxedSecondDerivative,
    const bool IsRelaxedSecondDerivativeHistorical)
    : mBossakAlpha(BossakAlpha),
      mOneMinusBossakAlpha(1.0 - BossakAlpha),
      mrSecondDerivative(rSecondDerivative),
      mrRelaxedSecondDerivative(rRelaxedSecondDerivative),
      mIsRelaxedSecondDerivativeHistorical(IsRelaxedSecondDerivativeHistorical)
{
    if (mIsRelaxedSecondDerivativeHistorical) {
        mRelaxedSecondDerivativeGetterMethod = &BossakPreviousStepSecondDerivativeCalculator::GetRelaxedSecondDerivativeFromHistorical;
    } else {
        mRelaxedSecondDerivativeGetterMethod = &BossakPreviousStepSecondDerivativeCalculator::GetRelaxedSecondDerivativeFromNonHistorical;
    }
}

template<class TVariableDataType>
void RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::Check(
    const ModelPart& rModelPart) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(mrSecondDerivative))
        << "Second derivative variable " << mrSecondDerivative.Name()
        << " not found in nodal solution step variables list in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF(mIsRelaxedSecondDerivativeHistorical &&  !rModelPart.HasNodalSolutionStepVariable(mrRelaxedSecondDerivative))
        << "Relaxed second derivative variable " << mrRelaxedSecondDerivative.Name()
        << " not found in nodal solution step variables list in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF(mBossakAlpha == 0.0)
        << "Bossak alpha is zero. Please provide a non-zero bossak alpha.\n";

    KRATOS_CATCH("");
}

template<class TVariableDataType>
void RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::Calculate(
    NodeType& rNode) const
{
    // get the relaxed second derivative
    const TVariableDataType& relaxed_second_derivative = (this->*(this->mRelaxedSecondDerivativeGetterMethod))(rNode);

    // compute the previous time step second derivative
    const TVariableDataType& r_current_second_derivative = rNode.FastGetSolutionStepValue(mrSecondDerivative, 0);
    rNode.FastGetSolutionStepValue(mrSecondDerivative, 1) = (relaxed_second_derivative - mOneMinusBossakAlpha * r_current_second_derivative) / mBossakAlpha;
}

template<class TVariableDataType>
const TVariableDataType& RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::GetRelaxedSecondDerivativeFromHistorical(
    const NodeType& rNode) const
{
    return rNode.FastGetSolutionStepValue(mrRelaxedSecondDerivative);
}

template<class TVariableDataType>
const TVariableDataType& RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::GetRelaxedSecondDerivativeFromNonHistorical(
    const NodeType& rNode) const
{
    return rNode.GetValue(mrRelaxedSecondDerivative);
}

template<class TVariableDataType>
std::string RansComputePreviousStepBossakSecondDerivativesProcess::BossakPreviousStepSecondDerivativeCalculator<TVariableDataType>::Info() const
{
    const auto& type_str = [](const bool IsHistorical) -> std::string {
        return (IsHistorical ? "Historical" : "Non-historical");
    };

    return "Historical " + mrSecondDerivative.Name() + " and " +
           type_str(mIsRelaxedSecondDerivativeHistorical) + " " +
           mrRelaxedSecondDerivative.Name();
}

RansComputePreviousStepBossakSecondDerivativesProcess::RansComputePreviousStepBossakSecondDerivativesProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    const double bossak_alpha = rParameters["bossak_alpha"].GetDouble();

    UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    const auto default_variable_data_paramters = Parameters(R"(
    {
        "second_derivative_variable_name"        : "PLEASE_SPECIFY_VARIABLE_NAME",
        "relaxed_second_derivative_variable_name": "PLEASE_SPECIFY_VARIABLE_NAME",
        "is_relaxed_second_derivative_historical": true
    })");

    auto variables_data_list_settings = rParameters["variables_list"];
    for (auto variable_data_settings : variables_data_list_settings) {
        variable_data_settings.ValidateAndAssignDefaults(default_variable_data_paramters);

        const auto& second_derivative_variable_name = variable_data_settings["second_derivative_variable_name"].GetString();
        const auto& relaxed_second_derivative_variable_name = variable_data_settings["relaxed_second_derivative_variable_name"].GetString();

        if (KratosComponents<Variable<double>>::Has(second_derivative_variable_name)) {
            if (KratosComponents<Variable<double>>::Has(relaxed_second_derivative_variable_name)) {
                mDoubleVariableDataList.push_back(BossakPreviousStepSecondDerivativeCalculator<double>(
                    bossak_alpha, KratosComponents<Variable<double>>::Get(second_derivative_variable_name),
                    KratosComponents<Variable<double>>::Get(relaxed_second_derivative_variable_name),
                    variable_data_settings["is_relaxed_second_derivative_historical"].GetBool()));
            } else {
                KRATOS_ERROR
                    << "Second derivative variable "
                    << second_derivative_variable_name << " is of type double, and relaxed second derivative variable "
                    << relaxed_second_derivative_variable_name << " is not of type double. Both of these variables should be of same type.\n";
            }
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(second_derivative_variable_name)) {
            if (KratosComponents<Variable<array_1d<double, 3>>>::Has(relaxed_second_derivative_variable_name)) {
                mArray3DVariableDataList.push_back(BossakPreviousStepSecondDerivativeCalculator<array_1d<double, 3>>(
                    bossak_alpha, KratosComponents<Variable<array_1d<double, 3>>>::Get(second_derivative_variable_name),
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(relaxed_second_derivative_variable_name),
                    variable_data_settings["is_relaxed_second_derivative_historical"].GetBool()));
            } else {
                KRATOS_ERROR
                    << "Second derivative variable "
                    << second_derivative_variable_name << " is of type array_1d<double, 3>, and relaxed second derivative variable "
                    << relaxed_second_derivative_variable_name << " is not of type array_1d<double, 3>. Both of these variables should be of same type.\n";
            }
        } else {
            KRATOS_ERROR
                << "Supports only double and array_1d<double, 3> variable "
                   "types. [ Second derivative varable name = "
                << second_derivative_variable_name << ", relaxed second derivatives variable name = "
                << relaxed_second_derivative_variable_name << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

int RansComputePreviousStepBossakSecondDerivativesProcess::Check()
{
    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(r_model_part.GetBufferSize() < 2)
        << mModelPartName << " needs to have a buffer size equal or greater than 2. [ Buffer size = "
        << r_model_part.GetBufferSize() << " ].\n";

    for (const auto& r_variable_data : mDoubleVariableDataList) {
        r_variable_data.Check(r_model_part);
    }

    for (const auto& r_variable_data : mArray3DVariableDataList) {
        r_variable_data.Check(r_model_part);
    }

    return 0;
}

void RansComputePreviousStepBossakSecondDerivativesProcess::ExecuteInitialize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE) {
            ExecuteCalculation();
            break;
        }
    }
}

void RansComputePreviousStepBossakSecondDerivativesProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE_SOLUTION_STEP) {
            ExecuteCalculation();
            break;
        }
    }
}

void RansComputePreviousStepBossakSecondDerivativesProcess::Execute()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::EXECUTE) {
            ExecuteCalculation();
            break;
        }
    }
}

void RansComputePreviousStepBossakSecondDerivativesProcess::ExecuteFinalizeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE_SOLUTION_STEP) {
            ExecuteCalculation();
            break;
        }
    }
}

void RansComputePreviousStepBossakSecondDerivativesProcess::ExecuteFinalize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE) {
            ExecuteCalculation();
            break;
        }
    }
}

const Parameters RansComputePreviousStepBossakSecondDerivativesProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_SOURCE_MODEL_PART_NAME",
            "bossak_alpha"    : -0.3,
            "execution_points": ["initialize"],
            "echo_level"      : 0,
            "variables_list"  : [
                {
                    "second_derivative_variable_name"        : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "relaxed_second_derivative_variable_name": "PLEASE_SPECIFY_VARIABLE_NAME",
                    "is_relaxed_second_derivative_historical": true
                }
            ]
        })");

    return default_parameters;
}

std::string RansComputePreviousStepBossakSecondDerivativesProcess::Info() const
{
    return std::string("RansComputePreviousStepBossakSecondDerivativesProcess");
}

void RansComputePreviousStepBossakSecondDerivativesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansComputePreviousStepBossakSecondDerivativesProcess::PrintData(std::ostream& rOStream) const
{
}

void RansComputePreviousStepBossakSecondDerivativesProcess::ExecuteCalculation()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode) {
        for (const auto& variable_data : mDoubleVariableDataList) {
            variable_data.Calculate(rNode);
        }
        for (const auto& variable_data : mArray3DVariableDataList) {
            variable_data.Calculate(rNode);
        }
    });

    if (mEchoLevel > 0) {
        std::stringstream msg;
        msg << "Previous step bossak second derivatives was calculated in " << mModelPartName << " using:\n";

        for (const auto& variable_data : mDoubleVariableDataList) {
            msg << "\t" << variable_data.Info() << "\n";
        }

        for (const auto& variable_data : mArray3DVariableDataList) {
            msg << "\t" << variable_data.Info() << "\n";
        }

        KRATOS_INFO(this->Info()) << msg.str();
    }

    KRATOS_CATCH("");
}

void RansComputePreviousStepBossakSecondDerivativesProcess::UpdateExecutionPointsList(const std::vector<std::string>& rExecutionPointsList)
{
    KRATOS_TRY

    mExecutionPointsList.clear();

    for (const auto& execution_point : rExecutionPointsList) {
        if (execution_point == "initialize") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE);
        } else if (execution_point == "initialize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE_SOLUTION_STEP);
        } else if (execution_point == "execute") {
            mExecutionPointsList.push_back(ExecutionPoint::EXECUTE);
        } else if (execution_point == "finalize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE_SOLUTION_STEP);
        } else if (execution_point == "finalize") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE);
        } else {
            KRATOS_ERROR << "Unsupported execution point provided. [ execution_point = " << execution_point << " ]. Supported points are: \n\n"
                        << "\tinitialize\n"
                        << "\tinitialize_solution_step\n"
                        << "\texecute\n"
                        << "\tfinalize_solution_step\n"
                        << "\tfinalize\n";
        }
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
