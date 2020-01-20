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
//                   Riccardo Tosi
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
            "averaged_variables_list"                       : [],
            "time_averaging_container"                      : "NodalHistorical",
            "time_averaging_method"                         : "Average",
            "integration_start_point_control_variable_name" : "TIME",
            "integration_start_point_control_value"         : 0.0,
            "echo_level"                                    : 0
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mVariableNamesList = mrParameters["variables_list"].GetStringArray();
    mAveragedVariableNamesList = mrParameters["averaged_variables_list"].GetStringArray();

    if (mrParameters["time_averaging_container"].GetString() ==
        "NodalHistorical")
        mTimeAveragingContainer = NodalHistorical;
    else if (mrParameters["time_averaging_container"].GetString() ==
             "NodalNonHistorical")
        mTimeAveragingContainer = NodalNonHistorical;
    else if (mrParameters["time_averaging_container"].GetString() ==
             "ElementalNonHistorical")
        mTimeAveragingContainer = ElementalNonHistorical;

    if (mrParameters["time_averaging_method"].GetString() == "Average")
        mTimeAveragingMethod = Average;
    else if (mrParameters["time_averaging_method"].GetString() ==
             "RootMeanSquare")
        mTimeAveragingMethod = RootMeanSquare;

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

    for (const std::string& variable_name : mVariableNamesList)
    {
        if (!(KratosComponents<Variable<double>>::Has(variable_name) ||
              KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)))
        {
            KRATOS_ERROR << "Variable " << variable_name << " not found in the scalar and 3d variables list of Kratos.\n";
        }
    }

    for (const std::string& variable_name : mAveragedVariableNamesList)
    {
        if (!(KratosComponents<Variable<double>>::Has(variable_name) ||
              KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)))
        {
            KRATOS_ERROR << "Variable " << variable_name << " not found in the scalar and 3d variables list of Kratos.\n";
        }
    }

    if (!(mAveragedVariableNamesList.size() == mVariableNamesList.size()))
    {
        KRATOS_ERROR << "List of input variables" << mVariableNamesList
                     << " and list of averaged variables"
                     << mAveragedVariableNamesList << " has different length. Please provide arrays of same length and corresponding variables of the same type.\n";
    }

    for (int counter = 0; counter < static_cast<int>(mVariableNamesList.size()); ++counter)
    {
        const std::string& variable_name = mVariableNamesList[counter];
        const std::string& averaged_variable_name = mAveragedVariableNamesList[counter];

        if (!((KratosComponents<Variable<double>>::Has(variable_name) &&
               KratosComponents<Variable<double>>::Has(averaged_variable_name)) ||
              (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name) &&
               KratosComponents<Variable<array_1d<double, 3>>>::Has(averaged_variable_name))))
        {
            KRATOS_ERROR << variable_name << " type and its corresponding averaging variable "
                         << averaged_variable_name << " has different types. Please provide corresponding variables of the same type.\n";
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
    ModelPart::ElementsContainerType& r_elements = r_communicator.LocalMesh().Elements();

    mCurrentTime = 0.0;

    std::stringstream msg;

    if (mTimeAveragingContainer == NodalHistorical || mTimeAveragingContainer == NodalNonHistorical)
    {
        for (const std::string& variable_name : mAveragedVariableNamesList)
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
    }
    else if (mTimeAveragingContainer == ElementalNonHistorical)
    {
        for (const std::string& variable_name : mAveragedVariableNamesList)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                VariableUtils().SetNonHistoricalVariableToZero(r_variable, r_elements);
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
                VariableUtils().SetNonHistoricalVariableToZero(r_variable, r_elements);
            }
            msg << " " << variable_name;
        }
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
        ModelPart::ElementsContainerType& r_elements =
            r_communicator.LocalMesh().Elements();

        const double delta_time = r_model_part.GetProcessInfo()[DELTA_TIME];

        std::stringstream msg;
        msg << "Integrating";

        for (int counter = 0; counter < static_cast<int>(mVariableNamesList.size()); ++counter)
        {
            const std::string& variable_name = mVariableNamesList[counter];
            const std::string& averaged_variable_name =
                mAveragedVariableNamesList[counter];

            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                const Variable<double>& r_averaged_variable =
                    KratosComponents<Variable<double>>::Get(averaged_variable_name);
                if (mTimeAveragingContainer == NodalHistorical)
                {
                    this->CalculateTimeIntegratedNodalHistoricalQuantity(
                        r_nodes, r_variable, r_averaged_variable, delta_time);
                }
                else if (mTimeAveragingContainer == NodalNonHistorical)
                {
                    this->CalculateTimeIntegratedNonHistoricalQuantity(
                        r_nodes, r_variable, r_averaged_variable, delta_time);
                }
                else if (mTimeAveragingContainer == ElementalNonHistorical)
                {
                    this->CalculateTimeIntegratedNonHistoricalQuantity(
                        r_elements, r_variable, r_averaged_variable, delta_time);
                }
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
                const Variable<array_1d<double, 3>>& r_averaged_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(averaged_variable_name);
                if (mTimeAveragingContainer == NodalHistorical)
                {
                    this->CalculateTimeIntegratedNodalHistoricalQuantity(
                        r_nodes, r_variable, r_averaged_variable, delta_time);
                }
                else if (mTimeAveragingContainer == NodalNonHistorical)
                {
                    this->CalculateTimeIntegratedNonHistoricalQuantity(
                        r_nodes, r_variable, r_averaged_variable, delta_time);
                }
                else if (mTimeAveragingContainer == ElementalNonHistorical)
                {
                    this->CalculateTimeIntegratedNonHistoricalQuantity(
                        r_elements, r_variable, r_averaged_variable, delta_time);
                }
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
std::function<void(const TDataType&, TDataType&, const double)> TimeAveragingProcess::GetTimeAveragingMethod() const
{
    if (mTimeAveragingMethod == Average)
    {
        return [this](const TDataType& rTimeValue, TDataType& rAveragedValue,
                      const double DeltaTime) {
            this->AverageMethod(rTimeValue, rAveragedValue, DeltaTime);
        };
    }
    else if (mTimeAveragingMethod == RootMeanSquare)
    {
        return [this](const TDataType& rTimeValue, TDataType& rAveragedValue,
                      const double DeltaTime) {
            this->RootMeanSquareMethod(rTimeValue, rAveragedValue, DeltaTime);
        };
    }
    else
    {
        KRATOS_ERROR << "Unsupported time averaging method.";
        return [this](const TDataType& rTimeValue, TDataType& rAveragedValue,
                      const double DeltaTime) {
            this->RootMeanSquareMethod(rTimeValue, rAveragedValue, DeltaTime);
        };
    }
}

template <typename TDataType>
void TimeAveragingProcess::CalculateTimeIntegratedNodalHistoricalQuantity(
    ModelPart::NodesContainerType& rNodes,
    const Variable<TDataType>& rVariable,
    const Variable<TDataType>& rAveragedVariable,
    const double DeltaTime) const
{
    const std::function<void(const TDataType&, TDataType&, const double)> averaging_method =
        this->GetTimeAveragingMethod<TDataType>();

    const int number_of_nodes = rNodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rNodes.begin() + i_node);
        const TDataType& r_temporal_value = r_node.FastGetSolutionStepValue(rVariable);
        TDataType& r_integrated_value = r_node.GetValue(rAveragedVariable);
        averaging_method(r_temporal_value, r_integrated_value, DeltaTime);
    }
}

template <typename TDataType, typename TContainerType>
void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity(
    TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const Variable<TDataType>& rAveragedVariable,
    const double DeltaTime) const
{
    const std::function<void(const TDataType&, TDataType&, const double)> averaging_method =
        this->GetTimeAveragingMethod<TDataType>();

    const int number_of_nodes = rContainer.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        auto& r_container_element = *(rContainer.begin() + i_node);
        const TDataType& r_temporal_value = r_container_element.GetValue(rVariable);
        TDataType& r_integrated_value = r_container_element.GetValue(rAveragedVariable);
        averaging_method(r_temporal_value, r_integrated_value, DeltaTime);
    }
}

template <typename TDataType>
void TimeAveragingProcess::AverageMethod(const TDataType& rTemporalVariable,
                                         TDataType& rAveragedVariable,
                                         const double DeltaTime) const
{
    rAveragedVariable = (rAveragedVariable * mCurrentTime + rTemporalVariable * DeltaTime) /
                        (mCurrentTime + DeltaTime);
}

template <>
void TimeAveragingProcess::RootMeanSquareMethod<array_1d<double, 3>>(
    const array_1d<double, 3>& rTemporalVariable,
    array_1d<double, 3>& rAveragedVariable,
    const double DeltaTime) const
{
    for (int i = 0; i < 3; ++i)
    {
        rAveragedVariable[i] =
            std::sqrt((std::pow(rAveragedVariable[i], 2) * mCurrentTime +
                       std::pow(rTemporalVariable[i], 2) * DeltaTime) /
                      (mCurrentTime + DeltaTime));
    }
}

template <typename TDataType>
void TimeAveragingProcess::RootMeanSquareMethod(const TDataType& rTemporalVariable,
                                                TDataType& rAveragedVariable,
                                                const double DeltaTime) const
{
    rAveragedVariable = std::sqrt((std::pow(rAveragedVariable, 2) * mCurrentTime +
                                   std::pow(rTemporalVariable, 2) * DeltaTime) /
                                  (mCurrentTime + DeltaTime));
}

// template instantiations

template void TimeAveragingProcess::CalculateTimeIntegratedNodalHistoricalQuantity<double>(
    ModelPart::NodesContainerType&, const Variable<double>&, const Variable<double>&, const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedNodalHistoricalQuantity<array_1d<double, 3>>(
    ModelPart::NodesContainerType&,
    const Variable<array_1d<double, 3>>&,
    const Variable<array_1d<double, 3>>&,
    const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity<double, ModelPart::NodesContainerType>(
    ModelPart::NodesContainerType&, const Variable<double>&, const Variable<double>&, const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity<array_1d<double, 3>, ModelPart::NodesContainerType>(
    ModelPart::NodesContainerType&,
    const Variable<array_1d<double, 3>>&,
    const Variable<array_1d<double, 3>>&,
    const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity<double, ModelPart::ElementsContainerType>(
    ModelPart::ElementsContainerType&,
    const Variable<double>&,
    const Variable<double>&,
    const double) const;

template void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity<array_1d<double, 3>, ModelPart::ElementsContainerType>(
    ModelPart::ElementsContainerType&,
    const Variable<array_1d<double, 3>>&,
    const Variable<array_1d<double, 3>>&,
    const double) const;

} // namespace Kratos.
