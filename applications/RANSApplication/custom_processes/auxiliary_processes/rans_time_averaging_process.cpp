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
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variables_list"  : [],
            "echo_level"      : 0
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mVariableNamesList = mrParameters["variables_list"].GetStringArray();

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
    }

    KRATOS_CATCH("");
}

void RansTimeAveragingProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    Communicator& r_communicator = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    const double delta_time = r_model_part.GetProcessInfo()[DELTA_TIME];
    mCurrentTime += delta_time;

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
    }

    KRATOS_CATCH("");
}

void RansTimeAveragingProcess::ExecuteFinalizeSolutionStep()
{
    Execute();
}

void RansTimeAveragingProcess::ExecuteFinalize()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    Communicator& r_communicator = r_model_part.GetCommunicator();
    ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    for (const std::string& variable_name : mVariableNamesList)
    {
        if (KratosComponents<Variable<double>>::Has(variable_name))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(variable_name);
            this->CalculateTimeAveragedQuantity(r_nodes, r_variable);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
            this->CalculateTimeAveragedQuantity(r_nodes, r_variable);
        }
    }

    KRATOS_CATCH("");
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
    ModelPart::NodesContainerType& rNodes, const Variable<TDataType>& rVariable)
{
    int number_of_nodes = rNodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        r_node.SetValue(rVariable, rVariable.Zero());
    }
}

template <typename TDataType>
void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity(
    ModelPart::NodesContainerType& rNodes, const Variable<TDataType>& rVariable, const double DeltaTime)
{
    int number_of_nodes = rNodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        const TDataType& r_temporal_value = r_node.FastGetSolutionStepValue(rVariable);
        TDataType& r_integrated_value = r_node.GetValue(rVariable);
        r_integrated_value += r_temporal_value * DeltaTime;
    }
}

template <typename TDataType>
void RansTimeAveragingProcess::CalculateTimeAveragedQuantity(
    ModelPart::NodesContainerType& rNodes, const Variable<TDataType>& rVariable)
{
    int number_of_nodes = rNodes.size();
    const double weighting_factor = 1.0 / mCurrentTime;

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        TDataType& r_integrated_value = r_node.GetValue(rVariable);
        r_integrated_value = r_integrated_value * weighting_factor;
    }
}

// template instantiations

template void RansTimeAveragingProcess::InitializeTimeAveragedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&);

template void RansTimeAveragingProcess::InitializeTimeAveragedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&);

template void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&, const double);

template void RansTimeAveragingProcess::CalculateTimeIntegratedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const double);

template void RansTimeAveragingProcess::CalculateTimeAveragedQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&);

template void RansTimeAveragingProcess::CalculateTimeAveragedQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&);

} // namespace Kratos.
