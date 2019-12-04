//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_time_averaging_process.h"

namespace Kratos
{
RansTimeAveragingProcess::RansTimeAveragingProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                               : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variables_list"                                : [],
            "integration_start_point_control_variable_name" : "TIME",
            "integration_start_point_control_value"         : 0.0,
            "echo_level"                                    : 0
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mVariableNamesList = mrParameters["variables_list"].GetStringArray();
    mIntegrationControlVariableName =
        mrParameters["integration_start_point_control_variable_name"].GetString();

    KRATOS_ERROR_IF(!KratosComponents<Variable<int>>::Has(mIntegrationControlVariableName) &&
                    !KratosComponents<Variable<double>>::Has(mIntegrationControlVariableName))
        << "\"integration_start_point_control_variable_name\" needs to be "
           "either integer or double variable name. Current variable name is \""
        << mIntegrationControlVariableName << "\".\n";

    KRATOS_CATCH("");
}

/// Destructor.
RansTimeAveragingProcess::~RansTimeAveragingProcess()
{
}

int RansTimeAveragingProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    for (const std::string& variable_name : mVariableNamesList)
    {
        if (KratosComponents<Variable<double>>::Has(variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(variable_name);
            RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, r_variable);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
            RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, r_variable);
        }
        else
        {
            KRATOS_ERROR << "Variable " << variable_name << " not found in the scalar and 3d variables list of Kratos.\n";
        }
    }

    return 0;

    KRATOS_CATCH("");
}

void RansTimeAveragingProcess::ExecuteInitialize()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    Communicator& r_communicator = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    mCurrentTime = 0.0;

    std::stringstream msg;

    msg << "Initialized non-historical";

    for (const std::string& variable_name : mVariableNamesList)
    {
        if (KratosComponents<Variable<double>>::Has(variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(variable_name);
            this->InitializeTimeAveragedQuantity(r_nodes, r_variable);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
            this->InitializeTimeAveragedQuantity(r_nodes, r_variable);
        }
        msg << " " << variable_name << ",";
    }

    msg.seekp(-1, msg.cur);
    msg << " variable(s) in " << mModelPartName << " to store time averaged quantities.\n";

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << msg.str();

    KRATOS_CATCH("");
}

void RansTimeAveragingProcess::Execute()
{
    KRATOS_TRY

    if (IsIntegrationStep())
    {
        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
        Communicator& r_communicator = r_model_part.GetCommunicator();
        ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

        const double delta_time = r_model_part.GetProcessInfo()[DELTA_TIME];

        std::stringstream msg;
        msg << "Integrating historical";

        for (const std::string& variable_name : mVariableNamesList)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                this->CalculateTimeIntegratedQuantity(r_nodes, r_variable, delta_time);
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
                this->CalculateTimeIntegratedQuantity(r_nodes, r_variable, delta_time);
            }

            msg << " " << variable_name << ",";
        }

        mCurrentTime += delta_time;

        msg.seekp(-1, msg.cur);
        msg << " variable(s) in " << mModelPartName << ".\n";

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1) << msg.str();
    }

    KRATOS_CATCH("");
}

void RansTimeAveragingProcess::ExecuteFinalizeSolutionStep()
{
    Execute();
}

bool RansTimeAveragingProcess::IsIntegrationStep() const
{
    const ProcessInfo& r_process_info =
        mrModel.GetModelPart(mModelPartName).GetProcessInfo();

    if (KratosComponents<Variable<int>>::Has(mIntegrationControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mIntegrationControlVariableName);
        return (r_process_info[r_variable] >=
                mrParameters["integration_start_point_control_value"].GetInt());
    }
    else if (KratosComponents<Variable<double>>::Has(mIntegrationControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mIntegrationControlVariableName);
        return (r_process_info[r_variable] >=
                mrParameters["integration_start_point_control_value"].GetDouble());
    }
    else
    {
        return false;
    }
}

std::string RansTimeAveragingProcess::Info() const
{
    return std::string("RansTimeAveragingProcess");
}

void RansTimeAveragingProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansTimeAveragingProcess::PrintData(std::ostream& rOStream) const
{
}

template <typename TDataType>
void RansTimeAveragingProcess::InitializeTimeAveragedQuantity(
    ModelPart::NodesContainerType& rNodes, const Variable<TDataType>& rVariable) const
{
    const int number_of_nodes = rNodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        r_node.SetValue(rVariable, rVariable.Zero());
    }
}

template <typename TDataType>
void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity(ModelPart::NodesContainerType& rNodes,
                                                               const Variable<TDataType>& rVariable,
                                                               const double DeltaTime) const
{
    const int number_of_nodes = rNodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        const TDataType& r_temporal_value = r_node.FastGetSolutionStepValue(rVariable);
        TDataType& r_integrated_value = r_node.GetValue(rVariable);
        r_integrated_value =
            (r_integrated_value * mCurrentTime + r_temporal_value * DeltaTime) /
            (mCurrentTime + DeltaTime);
    }
}

// template instantiations

template void RansTimeAveragingProcess::InitializeTimeAveragedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&) const;

template void RansTimeAveragingProcess::InitializeTimeAveragedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&) const;

template void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&, const double) const;

template void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const double) const;

} // namespace Kratos.
