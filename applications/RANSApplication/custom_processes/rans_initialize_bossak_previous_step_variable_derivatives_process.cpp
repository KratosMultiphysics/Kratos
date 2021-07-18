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
#include "utilities/time_discretization.h"

// Include base h
#include "rans_initialize_bossak_previous_step_variable_derivatives_process.h"

namespace Kratos
{
template<class TVariableDataType>
RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::BossakPreviousStepDerivativeCalculator(
    const Variable<TVariableDataType>& rFirstDerivativeVariable,
    const Variable<TVariableDataType>& rSecondDerivativeVariable,
    const Variable<TVariableDataType>& rRelaxedSecondDerivativeVariable,
    const bool IsRelaxedSecondDerivativeHistorical)
    : mrFirstDerivativeVariable(rFirstDerivativeVariable),
      mrSecondDerivativeVariable(rSecondDerivativeVariable),
      mrRelaxedSecondDerivativeVariable(rRelaxedSecondDerivativeVariable),
      mIsRelaxedSecondDerivativeVariableHistorical(IsRelaxedSecondDerivativeHistorical)
{
    if (mIsRelaxedSecondDerivativeVariableHistorical) {
        mRelaxedSecondDerivativeGetterMethod = &BossakPreviousStepDerivativeCalculator::GetRelaxedSecondDerivativeFromHistorical;
    } else {
        mRelaxedSecondDerivativeGetterMethod = &BossakPreviousStepDerivativeCalculator::GetRelaxedSecondDerivativeFromNonHistorical;
    }
}

template<class TVariableDataType>
void RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::Check(
    const ModelPart& rModelPart) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(mrFirstDerivativeVariable))
        << "First derivative variable " << mrFirstDerivativeVariable.Name()
        << " not found in nodal solution step variables list in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(mrSecondDerivativeVariable))
        << "Second derivative variable " << mrSecondDerivativeVariable.Name()
        << " not found in nodal solution step variables list in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF(mIsRelaxedSecondDerivativeVariableHistorical &&  !rModelPart.HasNodalSolutionStepVariable(mrRelaxedSecondDerivativeVariable))
        << "Relaxed second derivative variable " << mrRelaxedSecondDerivativeVariable.Name()
        << " not found in nodal solution step variables list in "
        << rModelPart.Name() << ".\n";

    KRATOS_CATCH("");
}

template<class TVariableDataType>
void RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::Calculate(
    NodeType& rNode,
    const BossakConstants& rBossakConstants) const
{
    // get the relaxed second derivative
    const TVariableDataType& relaxed_second_derivative = (this->*(this->mRelaxedSecondDerivativeGetterMethod))(rNode);

    const TVariableDataType& r_current_second_derivative = rNode.FastGetSolutionStepValue(mrSecondDerivativeVariable, 0);
    const TVariableDataType& r_current_first_derivative = rNode.FastGetSolutionStepValue(mrFirstDerivativeVariable, 0);

    TVariableDataType& r_old_second_derivative = rNode.FastGetSolutionStepValue(mrSecondDerivativeVariable, 1);
    r_old_second_derivative = (relaxed_second_derivative - rBossakConstants.C4 * r_current_second_derivative) / rBossakConstants.Alpha;

    TVariableDataType& r_old_first_derivative = rNode.FastGetSolutionStepValue(mrFirstDerivativeVariable, 1);
    r_old_first_derivative = r_current_first_derivative - (r_current_second_derivative + rBossakConstants.C3 * r_old_second_derivative) / rBossakConstants.C2;
}

template<class TVariableDataType>
const TVariableDataType& RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::GetRelaxedSecondDerivativeFromHistorical(
    const NodeType& rNode) const
{
    return rNode.FastGetSolutionStepValue(mrRelaxedSecondDerivativeVariable);
}

template<class TVariableDataType>
const TVariableDataType& RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::GetRelaxedSecondDerivativeFromNonHistorical(
    const NodeType& rNode) const
{
    return rNode.GetValue(mrRelaxedSecondDerivativeVariable);
}

template<class TVariableDataType>
std::string RansInitializeBossakPreviousStepVariableDerivatives::BossakPreviousStepDerivativeCalculator<TVariableDataType>::Info() const
{
    const auto& type_str = [](const bool IsHistorical) -> std::string {
        return (IsHistorical ? "Historical" : "Non-historical");
    };

    return "Historical " + mrFirstDerivativeVariable.Name() + ", Historical " + mrSecondDerivativeVariable.Name() + ", and " +
           type_str(mIsRelaxedSecondDerivativeVariableHistorical) + " " +
           mrRelaxedSecondDerivativeVariable.Name();
}

RansInitializeBossakPreviousStepVariableDerivatives::RansInitializeBossakPreviousStepVariableDerivatives(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mBossakAlpha = rParameters["bossak_alpha"].GetDouble();

    auto variables_data_list_settings = rParameters["variables_list"];
    for (auto variable_data_settings : variables_data_list_settings) {
        variable_data_settings.ValidateAndAssignDefaults(GetDefaultParameters()["variables_list"][0]);

        const auto& first_derivative_variable_name = variable_data_settings["first_derivative_variable_name"].GetString();
        const auto& second_derivative_variable_name = variable_data_settings["second_derivative_variable_name"].GetString();
        const auto& relaxed_second_derivative_variable_name = variable_data_settings["relaxed_second_derivative_variable_name"].GetString();

        if (CheckVariableTypes<double>(first_derivative_variable_name, second_derivative_variable_name, relaxed_second_derivative_variable_name)) {
            mDoubleVariableDataList.push_back(BossakPreviousStepDerivativeCalculator<double>(
                KratosComponents<Variable<double>>::Get(first_derivative_variable_name),
                KratosComponents<Variable<double>>::Get(second_derivative_variable_name),
                KratosComponents<Variable<double>>::Get(relaxed_second_derivative_variable_name),
                variable_data_settings["is_relaxed_second_derivative_historical"].GetBool()));
        } else if (CheckVariableTypes<array_1d<double, 3>>(first_derivative_variable_name, second_derivative_variable_name, relaxed_second_derivative_variable_name)) {
            mArray3DVariableDataList.push_back(BossakPreviousStepDerivativeCalculator<array_1d<double, 3>>(
                KratosComponents<Variable<array_1d<double, 3>>>::Get(first_derivative_variable_name),
                KratosComponents<Variable<array_1d<double, 3>>>::Get(second_derivative_variable_name),
                KratosComponents<Variable<array_1d<double, 3>>>::Get(relaxed_second_derivative_variable_name),
                variable_data_settings["is_relaxed_second_derivative_historical"].GetBool()));
        } else {
            KRATOS_ERROR
                << "Supports only double and array_1d<double, 3> variable types. [ "
                << "First derivative variable name = " << first_derivative_variable_name
                << ", Second derivative varable name = " << second_derivative_variable_name
                << ", relaxed second derivatives variable name = "
                << relaxed_second_derivative_variable_name << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

int RansInitializeBossakPreviousStepVariableDerivatives::Check()
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

void RansInitializeBossakPreviousStepVariableDerivatives::ExecuteInitializeSolutionStep()
{
    if (! mIsInitialized) {
        ExecuteCalculation();
        mIsInitialized = true;
    }
}

const Parameters RansInitializeBossakPreviousStepVariableDerivatives::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_SOURCE_MODEL_PART_NAME",
            "bossak_alpha"    : -0.3,
            "echo_level"      : 0,
            "variables_list"  : [
                {
                    "first_derivative_variable_name"         : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "second_derivative_variable_name"        : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "relaxed_second_derivative_variable_name": "PLEASE_SPECIFY_VARIABLE_NAME",
                    "is_relaxed_second_derivative_historical": true
                }
            ]
        })");

    return default_parameters;
}

std::string RansInitializeBossakPreviousStepVariableDerivatives::Info() const
{
    return std::string("RansInitializeBossakPreviousStepVariableDerivatives");
}

void RansInitializeBossakPreviousStepVariableDerivatives::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansInitializeBossakPreviousStepVariableDerivatives::PrintData(std::ostream& rOStream) const
{
}

void RansInitializeBossakPreviousStepVariableDerivatives::ExecuteCalculation()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const double delta_time = r_model_part.GetProcessInfo()[DELTA_TIME];

    TimeDiscretization::Bossak bossak(mBossakAlpha, 0.25, 0.5);

    BossakConstants bossak_constants;
    bossak_constants.Alpha = bossak.GetAlphaM();
    bossak_constants.Gamma = bossak.GetGamma();
    bossak_constants.C0 = (1.0 - bossak_constants.Alpha) / (bossak_constants.Gamma * delta_time);
    bossak_constants.C2 = 1.0 / (bossak_constants.Gamma * delta_time);
    bossak_constants.C3 = (1.0 - bossak_constants.Gamma) / bossak_constants.Gamma;
    bossak_constants.C4 = (1.0 - bossak_constants.Alpha);

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode) {
        for (const auto& variable_data : mDoubleVariableDataList) {
            variable_data.Calculate(rNode, bossak_constants);
        }
        for (const auto& variable_data : mArray3DVariableDataList) {
            variable_data.Calculate(rNode, bossak_constants);
        }
    });

    if (mEchoLevel > 0) {
        std::stringstream msg;
        msg << "Previous step bossak derivatives was calculated in " << mModelPartName << " using:\n";

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

} // namespace Kratos.
