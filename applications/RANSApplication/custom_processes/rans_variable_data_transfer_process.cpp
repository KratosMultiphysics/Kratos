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
#include "rans_variable_data_transfer_process.h"

namespace Kratos
{
template<class TVariableType>
RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyVariableData(
    const TVariableType& rSourceVariable,
    const bool IsSourceVariableHistorical,
    const IndexType SourceStepIndex,
    const TVariableType& rDestinationVariable,
    const bool IsDestinationVariableHistorical,
    const IndexType DestinationStepIndex)
    : mrSourceVariable(rSourceVariable),
      mIsSourceVariableHistorical(IsSourceVariableHistorical),
      mSourceStepIndex(SourceStepIndex),
      mrDestinationVariable(rDestinationVariable),
      mIsDestinationVariableHistorical(IsDestinationVariableHistorical),
      mDestinationStepIndex(DestinationStepIndex)
{
    if (mIsSourceVariableHistorical) {
        if (mIsDestinationVariableHistorical) {
            mCopyMethod = &CopyVariableData::CopyHistoricalToHistorical;
        } else {
            mCopyMethod = &CopyVariableData::CopyHistoricalToNonHistorical;
        }
    } else {
        if (mIsDestinationVariableHistorical) {
            mCopyMethod = &CopyVariableData::CopyNonHistoricalToHistorical;
        } else {
            mCopyMethod = &CopyVariableData::CopyNonHistoricalToNonHistorical;
        }
    }
}

template<class TVariableType>
void RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyData(
    const NodeType& rSourceNode,
    NodeType& rDestinationNode) const
{
    (this->*(this->mCopyMethod))(rSourceNode, rDestinationNode);
}

template<class TVariableType>
void RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.FastGetSolutionStepValue(mrDestinationVariable, mDestinationStepIndex) = rSourceNode.FastGetSolutionStepValue(mrSourceVariable, mSourceStepIndex);
}

template<class TVariableType>
void RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.SetValue(mrDestinationVariable, rSourceNode.FastGetSolutionStepValue(mrSourceVariable, mSourceStepIndex));
}

template<class TVariableType>
void RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyNonHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.FastGetSolutionStepValue(mrDestinationVariable, mDestinationStepIndex) = rSourceNode.GetValue(mrSourceVariable);
}

template<class TVariableType>
void RansVariableDataTransferProcess::CopyVariableData<TVariableType>::CopyNonHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.SetValue(mrDestinationVariable, rSourceNode.GetValue(mrSourceVariable));
}

template<class TVariableType>
std::string RansVariableDataTransferProcess::CopyVariableData<TVariableType>::Info() const
{
    const auto& type_str = [](const bool IsHistorical, const int StepIndex) -> std::string {
        return (IsHistorical ? ("Historical at step " + std::to_string(StepIndex)) : "Non-historical");
    };

    return type_str(mIsSourceVariableHistorical, mSourceStepIndex) + " " + mrSourceVariable.Name() +
           " to " + type_str(mIsDestinationVariableHistorical, mDestinationStepIndex) + " " +
           mrDestinationVariable.Name();
}

template<class TVariableType>
void RansVariableDataTransferProcess::CheckVariableData(
    const ModelPart& rModelPart,
    const bool CheckSourceVariableData,
    const CopyVariableData<TVariableType>& rCopyVariableData) const
{
    KRATOS_TRY

    if (CheckSourceVariableData) {
        if (rCopyVariableData.mIsSourceVariableHistorical) {
            KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(
                rCopyVariableData.mrSourceVariable))
                << "Source " << rCopyVariableData.mrSourceVariable.Name()
                << " is not found in solution step variables list of "
                << rModelPart.FullName() << ".\n";
            KRATOS_ERROR_IF(rModelPart.GetBufferSize() <= rCopyVariableData.mSourceStepIndex)
                << "Requested source step index of " << rCopyVariableData.mSourceStepIndex
                << " for variable " << rCopyVariableData.mrSourceVariable.Name() << " is not available in "
                << rModelPart.FullName() << ". Maximum allowed step index is "
                << rModelPart.GetBufferSize() - 1 << ".\n";
        } else {
            KRATOS_ERROR_IF(rCopyVariableData.mSourceStepIndex != 0)
                << "Source step index should be zero for non-historical variables. "
                "[ Source variable name = "
                << rCopyVariableData.mrSourceVariable.Name()
                << ", step index = " << rCopyVariableData.mSourceStepIndex << " ].\n";
        }
    } else {
        if (rCopyVariableData.mIsDestinationVariableHistorical) {
            KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(
                rCopyVariableData.mrDestinationVariable))
                << "Destination " << rCopyVariableData.mrDestinationVariable.Name()
                << " is not found in solution step variables list of "
                << rModelPart.FullName() << ".\n";
            KRATOS_ERROR_IF(rModelPart.GetBufferSize() <= rCopyVariableData.mDestinationStepIndex)
                << "Requested destination step index of " << rCopyVariableData.mDestinationStepIndex
                << " for variable " << rCopyVariableData.mrDestinationVariable.Name() << " is not available in "
                << rModelPart.FullName() << ". Maximum allowed step index is "
                << rModelPart.GetBufferSize() - 1 << ".\n";
        } else {
            KRATOS_ERROR_IF(rCopyVariableData.mDestinationStepIndex != 0)
                << "Destination step index should be zero for non-historical variables. "
                "[ Destination variable name = "
                << rCopyVariableData.mrDestinationVariable.Name()
                << ", step index = " << rCopyVariableData.mDestinationStepIndex << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

RansVariableDataTransferProcess::RansVariableDataTransferProcess(
    Model& rModel,
    Parameters rParameters)
    : mrSourceModel(rModel),
      mrDestinationModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mSourceModelPartName = rParameters["source_model_part_name"].GetString();
    mDestinationModelPartName = rParameters["destination_model_part_name"].GetString();

    const auto& copy_execution_points = rParameters["copy_execution_points"].GetStringArray();
    this->UpdateExecutionPointsList(copy_execution_points);

    const auto default_copy_variable_data_paramters = Parameters(R"(
    {
        "source_variable_settings" : {
            "is_historical_container": true,
            "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME",
            "step_index"             : 0
        },
        "destination_variable_settings" : {
            "is_historical_container": false,
            "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME",
            "step_index"             : 0
        }
    })");

    auto copy_variables_data_list_settings = rParameters["copy_variables_list"];
    std::vector<CopyVariableDataListItem> copy_variables_data_list;
    for (auto copy_variable_data_settings : copy_variables_data_list_settings) {
        copy_variable_data_settings.RecursivelyValidateAndAssignDefaults(default_copy_variable_data_paramters);

        const auto& source_variable_settings = copy_variable_data_settings["source_variable_settings"];
        const auto& destination_variable_settings = copy_variable_data_settings["destination_variable_settings"];

        copy_variables_data_list.push_back(std::make_tuple(
            source_variable_settings["variable_name"].GetString(),
            source_variable_settings["is_historical_container"].GetBool(),
            source_variable_settings["step_index"].GetInt(),
            destination_variable_settings["variable_name"].GetString(),
            destination_variable_settings["is_historical_container"].GetBool(),
            destination_variable_settings["step_index"].GetInt()
        ));
    }

    UpdateCopyVariableDataList(copy_variables_data_list);

    KRATOS_CATCH("");
}

RansVariableDataTransferProcess::RansVariableDataTransferProcess(
    Model& rModel,
    const std::string& rSourceModelPartName,
    const std::string& rDestinationModelPartName,
    const std::vector<std::string>& rCopyExecutionPoints,
    const std::vector<CopyVariableDataListItem>& rCopyVariableDataList,
    const int EchoLevel)
    : mrSourceModel(rModel),
      mrDestinationModel(rModel),
      mEchoLevel(EchoLevel),
      mSourceModelPartName(rSourceModelPartName),
      mDestinationModelPartName(rDestinationModelPartName)
{
    KRATOS_TRY

    this->UpdateExecutionPointsList(rCopyExecutionPoints);

    UpdateCopyVariableDataList(rCopyVariableDataList);

    KRATOS_CATCH("");
}

RansVariableDataTransferProcess::RansVariableDataTransferProcess(
    Model& rSourceModel,
    Model& rDestinationModel,
    const std::string& rSourceModelPartName,
    const std::string& rDestinationModelPartName,
    const std::vector<std::string>& rCopyExecutionPoints,
    const std::vector<CopyVariableDataListItem>& rCopyVariableDataList,
    const int EchoLevel)
    : mrSourceModel(rSourceModel),
      mrDestinationModel(rDestinationModel),
      mEchoLevel(EchoLevel),
      mSourceModelPartName(rSourceModelPartName),
      mDestinationModelPartName(rDestinationModelPartName)
{
    KRATOS_TRY

    this->UpdateExecutionPointsList(rCopyExecutionPoints);

    UpdateCopyVariableDataList(rCopyVariableDataList);

    KRATOS_CATCH("");
}

int RansVariableDataTransferProcess::Check()
{
    const auto& r_source_model_part = mrSourceModel.GetModelPart(mSourceModelPartName);
    const auto& r_source_model_part_nodes = r_source_model_part.Nodes();
    const auto& r_destination_model_part = mrDestinationModel.GetModelPart(mDestinationModelPartName);

    KRATOS_ERROR_IF(r_source_model_part_nodes.size() !=
                    r_destination_model_part.Nodes().size())
        << "Number of nodes mismatch in sourc " << mSourceModelPartName
        << " and destination " << mDestinationModelPartName
        << " [ number of nodes in source = " << r_source_model_part_nodes.size()
        << ", number of nodes in destination model part = "
        << r_destination_model_part.Nodes().size() << " ].\n";

    for (const auto& copy_variable_data : mCopyDoubleVariableDataList) {
        CheckVariableData(r_source_model_part, true, copy_variable_data);
        CheckVariableData(r_destination_model_part, false, copy_variable_data);
    }

    for (const auto& copy_variable_data : mCopyArray3DVariableDataList) {
        CheckVariableData(r_source_model_part, true, copy_variable_data);
        CheckVariableData(r_destination_model_part, false, copy_variable_data);
    }

    for (const auto& copy_variable_data : mCopyVectorVariableDataList) {
        CheckVariableData(r_source_model_part, true, copy_variable_data);
        CheckVariableData(r_destination_model_part, false, copy_variable_data);
    }

    for (const auto& copy_variable_data : mCopyMatrixVariableDataList) {
        CheckVariableData(r_source_model_part, true, copy_variable_data);
        CheckVariableData(r_destination_model_part, false, copy_variable_data);
    }

    return 0;
}

const Parameters RansVariableDataTransferProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "source_model_part_name"     : "PLEASE_SPECIFY_SOURCE_MODEL_PART_NAME",
            "destination_model_part_name": "PLEASE_SPECIFY_DESTINATION_MODEL_PART_NAME",
            "copy_execution_points"      : ["initialize"],
            "copy_variables_list"        : [
                {
                    "source_variable_settings" : {
                        "is_historical_container": true,
                        "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME",
                        "step_index"             : 0
                    },
                    "destination_variable_settings" : {
                        "is_historical_container": false,
                        "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME",
                        "step_index"             : 0
                    }
                }
            ],
            "echo_level"                 : 0
        })");

    return default_parameters;
}

std::string RansVariableDataTransferProcess::Info() const
{
    return std::string("RansVariableDataTransferProcess");
}

void RansVariableDataTransferProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansVariableDataTransferProcess::PrintData(std::ostream& rOStream) const
{
}

void RansVariableDataTransferProcess::ExecuteOperation()
{
    KRATOS_TRY

    const auto& r_source_model_part_nodes = mrSourceModel.GetModelPart(mSourceModelPartName).Nodes();
    auto& r_destination_model_part = mrDestinationModel.GetModelPart(mDestinationModelPartName);

    KRATOS_ERROR_IF(r_source_model_part_nodes.size() !=
                    r_destination_model_part.Nodes().size())
        << "Number of nodes mismatch in sourc " << mSourceModelPartName
        << " and destination " << mDestinationModelPartName
        << " [ number of nodes in source = " << r_source_model_part_nodes.size()
        << ", number of nodes in destination model part = "
        << r_destination_model_part.Nodes().size() << " ].\n";

    block_for_each(r_source_model_part_nodes, [&](const NodeType& rSourceNode) {
        auto& r_destination_node = r_destination_model_part.GetNode(rSourceNode.Id());
        for (const auto& copy_variable_data : mCopyDoubleVariableDataList) {
            copy_variable_data.CopyData(rSourceNode, r_destination_node);
        }
        for (const auto& copy_variable_data : mCopyArray3DVariableDataList) {
            copy_variable_data.CopyData(rSourceNode, r_destination_node);
        }
        for (const auto& copy_variable_data : mCopyVectorVariableDataList) {
            copy_variable_data.CopyData(rSourceNode, r_destination_node);
        }
        for (const auto& copy_variable_data : mCopyMatrixVariableDataList) {
            copy_variable_data.CopyData(rSourceNode, r_destination_node);
        }
    });

    // In here no mpi communication isrequired because we assume source model part is properly synchronized.

    if (mEchoLevel > 0) {
        std::stringstream msg;
        msg << "Following variable data was copied successfully from "
            << mSourceModelPartName << " to " << mDestinationModelPartName << ":\n";

        for (const auto& copy_variable_data : mCopyDoubleVariableDataList) {
            msg << "\t" << copy_variable_data.Info() << "\n";
        }

        for (const auto& copy_variable_data : mCopyArray3DVariableDataList) {
            msg << "\t" << copy_variable_data.Info() << "\n";
        }

        for (const auto& copy_variable_data : mCopyVectorVariableDataList) {
            msg << "\t" << copy_variable_data.Info() << "\n";
        }

        for (const auto& copy_variable_data : mCopyMatrixVariableDataList) {
            msg << "\t" << copy_variable_data.Info() << "\n";
        }

        KRATOS_INFO(this->Info()) << msg.str();
    }

    KRATOS_CATCH("");
}

void RansVariableDataTransferProcess::UpdateCopyVariableDataList(const std::vector<CopyVariableDataListItem>& rCopyVariableDataList)
{
    KRATOS_TRY

    mCopyDoubleVariableDataList.clear();
    mCopyArray3DVariableDataList.clear();
    mCopyVectorVariableDataList.clear();
    mCopyMatrixVariableDataList.clear();

    for (const auto&  r_copy_variable_data_item : rCopyVariableDataList) {
        const auto& source_variable_name = std::get<0>(r_copy_variable_data_item);
        const bool is_source_variable_historical = std::get<1>(r_copy_variable_data_item);
        const int source_step_index = std::get<2>(r_copy_variable_data_item);
        const auto& destination_variable_name = std::get<3>(r_copy_variable_data_item);
        const bool is_destination_variable_historical = std::get<4>(r_copy_variable_data_item);
        const int destination_step_index = std::get<5>(r_copy_variable_data_item);

        if (KratosComponents<Variable<double>>::Has(source_variable_name)) {
            const auto& r_source_variable = KratosComponents<Variable<double>>::Get(source_variable_name);

            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(destination_variable_name))
                << "Source variable " << source_variable_name
                << " is of double variable type, but destination variable "
                << destination_variable_name
                << " not found in double variables database." << std::endl;
            const auto& r_destination_variable = KratosComponents<Variable<double>>::Get(destination_variable_name);

            mCopyDoubleVariableDataList.push_back(CopyVariableData<Variable<double>>(
                r_source_variable, is_source_variable_historical, source_step_index,
                r_destination_variable, is_destination_variable_historical, destination_step_index));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(source_variable_name)) {
            const auto& r_source_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(source_variable_name);
            const bool has_destination_variable = KratosComponents<Variable<array_1d<double, 3>>>::Has(destination_variable_name);

            KRATOS_ERROR_IF_NOT(has_destination_variable)
                << "Source variable " << source_variable_name << " is of array_1d<double, 3> variable type, but destination variable "
                << destination_variable_name
                << " not found in array_1d<double, 3> variables database." << std::endl;
            const auto& r_destination_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(destination_variable_name);

            mCopyArray3DVariableDataList.push_back(CopyVariableData<Variable<array_1d<double, 3>>>(
                r_source_variable, is_source_variable_historical, source_step_index,
                r_destination_variable, is_destination_variable_historical, destination_step_index));
        } else if (KratosComponents<Variable<Vector>>::Has(source_variable_name)) {
            const auto& r_source_variable = KratosComponents<Variable<Vector>>::Get(source_variable_name);

            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<Vector>>::Has(destination_variable_name))
                << "Source variable " << source_variable_name
                << " is of Vector variable type, but destination variable "
                << destination_variable_name
                << " not found in Vector variables database." << std::endl;
            const auto& r_destination_variable = KratosComponents<Variable<Vector>>::Get(destination_variable_name);

            mCopyVectorVariableDataList.push_back(CopyVariableData<Variable<Vector>>(
                r_source_variable, is_source_variable_historical, source_step_index,
                r_destination_variable, is_destination_variable_historical, destination_step_index));
        } else if (KratosComponents<Variable<Matrix>>::Has(source_variable_name)) {
            const auto& r_source_variable = KratosComponents<Variable<Matrix>>::Get(source_variable_name);

            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<Matrix>>::Has(destination_variable_name))
                << "Source variable " << source_variable_name
                << " is of Matrix variable type, but destination variable "
                << destination_variable_name
                << " not found in Matrix variables database." << std::endl;
            const auto& r_destination_variable = KratosComponents<Variable<Matrix>>::Get(destination_variable_name);

            mCopyMatrixVariableDataList.push_back(CopyVariableData<Variable<Matrix>>(
                r_source_variable, is_source_variable_historical, source_step_index,
                r_destination_variable, is_destination_variable_historical, destination_step_index));
        } else {
            KRATOS_ERROR << "Unsupported source variable type requested. [ "
                            "variable name = "
                         << source_variable_name << " ]. Supported variable types are"
                         << "\n\t double"
                         << "\n\t array_1d<double, 3>"
                         << "\n\t Vector"
                         << "\n\t Matrix";
        }
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
