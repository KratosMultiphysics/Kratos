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
#include "utilities/variable_utils.h"

// Include base h
#include "time_averaging_process.h"

namespace Kratos
{
TimeAveragingProcess::TimeAveragingProcess(Model& rModel, Parameters rParameters)
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
TimeAveragingProcess::~TimeAveragingProcess()
{
}

int TimeAveragingProcess::Check()
{
    KRATOS_TRY

    if (!mrModel.HasModelPart(mModelPartName))
    {
        const std::vector<std::string>& r_model_part_names = mrModel.GetModelPartNames();

        std::string msg;
        msg = mrModel.Info() + " doesn't have " + mModelPartName +
              ". Available model parts are: \n";
        for (std::string model_part_name : r_model_part_names)
            msg += "     " + model_part_name + "\n";

        KRATOS_ERROR << msg;
    }

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    for (const std::string& variable_name : mVariableNamesList)
    {
        if (KratosComponents<Variable<double>>::Has(variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(variable_name);
            KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_variable))
                << r_model_part.Name() << " doesn't have "
                << r_variable.Name() << " in NodalSolutionStepDataContainer. Please add it as a SolutionStepVariable.";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
            KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_variable))
                << r_model_part.Name() << " doesn't have "
                << r_variable.Name() << " in NodalSolutionStepDataContainer. Please add it as a SolutionStepVariable.";
        }
        else
        {
            KRATOS_ERROR << "Variable " << variable_name << " not found in the scalar and 3d variables list of Kratos.\n";
        }
    }

    return 0;

    KRATOS_CATCH("");
}

void TimeAveragingProcess::ExecuteInitialize()
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
            VariableUtils().SetNonHistoricalVariableToZero(r_variable, r_nodes);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
            VariableUtils().SetNonHistoricalVariableToZero(r_variable, r_nodes);
        }
        msg << " " << variable_name;
    }
    msg << " variable(s) in " << mModelPartName << " to store time averaged quantities.\n";

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << msg.str();

    KRATOS_CATCH("");
}

void TimeAveragingProcess::ExecuteFinalizeSolutionStep()
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

            msg << " " << variable_name;
        }

        mCurrentTime += delta_time;

        msg << " variable(s) in " << mModelPartName << ".\n";
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1) << msg.str();
    }

    KRATOS_CATCH("");
}

bool TimeAveragingProcess::IsIntegrationStep() const
{
    const ProcessInfo& r_process_info =
        mrModel.GetModelPart(mModelPartName).GetProcessInfo();

    if (KratosComponents<Variable<int>>::Has(mIntegrationControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mIntegrationControlVariableName);

        KRATOS_ERROR_IF(!r_process_info.Has(r_variable))
            << "\"" << mIntegrationControlVariableName
            << "\" not found in process info of " << mModelPartName << ".\n";

        return (r_process_info[r_variable] >=
                mrParameters["integration_start_point_control_value"].GetInt());
    }
    else if (KratosComponents<Variable<double>>::Has(mIntegrationControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mIntegrationControlVariableName);

        KRATOS_ERROR_IF(!r_process_info.Has(r_variable))
            << "\"" << mIntegrationControlVariableName
            << "\" not found in process info of " << mModelPartName << ".\n";

        return (r_process_info[r_variable] >=
                mrParameters["integration_start_point_control_value"].GetDouble());
    }
    else
    {
        return false;
    }
}

std::string TimeAveragingProcess::Info() const
{
    return std::string("TimeAveragingProcess");
}

void TimeAveragingProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void TimeAveragingProcess::PrintData(std::ostream& rOStream) const
{
}

template <typename TDataType>
void TimeAveragingProcess::CalculateTimeIntegratedQuantity(ModelPart::NodesContainerType& rNodes,
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

template void TimeAveragingProcess::CalculateTimeIntegratedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&, const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const double) const;

} // namespace Kratos.