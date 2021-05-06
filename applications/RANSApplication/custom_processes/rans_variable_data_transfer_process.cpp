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
RansVariableDataTransferProcess::CopyVariableData::CopyVariableData(
    const VariableType& rSourceVariable,
    const bool IsSourceVariableHistorical,
    const VariableType& rDestinationVariable,
    const bool IsDestinationVariableHistorical)
    : mrSourceVariable(rSourceVariable),
      mIsSourceVariableHistorical(IsSourceVariableHistorical),
      mrDestinationVariable(rDestinationVariable),
      mIsDestinationVariableHistorical(IsDestinationVariableHistorical)
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

void RansVariableDataTransferProcess::CopyVariableData::CopyData(
    const NodeType& rSourceNode,
    NodeType& rDestinationNode) const
{
    (this->*(this->mCopyMethod))(rSourceNode, rDestinationNode);
}

void RansVariableDataTransferProcess::CopyVariableData::CopyHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.FastGetSolutionStepValue(mrDestinationVariable) = rSourceNode.FastGetSolutionStepValue(mrSourceVariable);
}

void RansVariableDataTransferProcess::CopyVariableData::CopyHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.SetValue(mrDestinationVariable, rSourceNode.FastGetSolutionStepValue(mrSourceVariable));
}

void RansVariableDataTransferProcess::CopyVariableData::CopyNonHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.FastGetSolutionStepValue(mrDestinationVariable) = rSourceNode.GetValue(mrSourceVariable);
}

void RansVariableDataTransferProcess::CopyVariableData::CopyNonHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const
{
    rDestinationNode.SetValue(mrDestinationVariable, rSourceNode.GetValue(mrSourceVariable));
}

std::string RansVariableDataTransferProcess::CopyVariableData::Info() const
{
    const auto& type_str = [](const bool IsHistorical) -> std::string {
        return (IsHistorical ? "Historical" : "Non-historical");
    };

    return type_str(mIsSourceVariableHistorical) + " " + mrSourceVariable.Name() +
           " to " + type_str(mIsDestinationVariableHistorical) + " " +
           mrDestinationVariable.Name();
}

RansVariableDataTransferProcess::RansVariableDataTransferProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mSourceModelPartName = rParameters["source_model_part_name"].GetString();
    mDestinationModelPartName = rParameters["destination_model_part_name"].GetString();

    const auto& copy_execution_points = rParameters["copy_execution_points"].GetStringArray();
    UpdateCopyExecutionPointsList(copy_execution_points);

    const auto default_copy_variable_data_paramters = Parameters(R"(
    {
        "source_variable_settings" : {
            "is_historical_container": true,
            "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME"
        },
        "destination_variable_settings" : {
            "is_historical_container": false,
            "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME"
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
            destination_variable_settings["variable_name"].GetString(),
            destination_variable_settings["is_historical_container"].GetBool()
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
    : mrModel(rModel),
      mEchoLevel(EchoLevel),
      mSourceModelPartName(rSourceModelPartName),
      mDestinationModelPartName(rDestinationModelPartName)
{
    KRATOS_TRY

    UpdateCopyExecutionPointsList(rCopyExecutionPoints);

    UpdateCopyVariableDataList(rCopyVariableDataList);

    KRATOS_CATCH("");
}

int RansVariableDataTransferProcess::Check()
{
    const auto& r_source_model_part_nodes = mrModel.GetModelPart(mSourceModelPartName).Nodes();
    auto& r_destination_model_part = mrModel.GetModelPart(mDestinationModelPartName);

    KRATOS_ERROR_IF(r_source_model_part_nodes.size() !=
                    r_destination_model_part.Nodes().size())
        << "Number of nodes mismatch in sourc " << mSourceModelPartName
        << " and destination " << mDestinationModelPartName
        << " [ number of nodes in source = " << r_source_model_part_nodes.size()
        << ", number of nodes in destination model part = "
        << r_destination_model_part.Nodes().size() << " ].\n";

    return 0;
}

void RansVariableDataTransferProcess::ExecuteInitialize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE_SOLUTION_STEP) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::ExecuteBeforeCouplingSolveStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::BEFORE_COUPLING_SOLVE_STEP) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::Execute()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::EXECUTE) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::ExecuteAfterCouplingSolveStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::AFTER_COUPLING_SOLVE_STEP) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::ExecuteFinalizeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE_SOLUTION_STEP) {
            ExecuteCopy();
            break;
        }
    }
}

void RansVariableDataTransferProcess::ExecuteFinalize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE) {
            ExecuteCopy();
            break;
        }
    }
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
                        "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME"
                    },
                    "destination_variable_settings" : {
                        "is_historical_container": false,
                        "variable_name"          : "PLEASE_SPECIFY_VARIABLE_NAME"
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

void RansVariableDataTransferProcess::ExecuteCopy()
{
    KRATOS_TRY

    const auto& r_source_model_part_nodes = mrModel.GetModelPart(mSourceModelPartName).Nodes();
    auto& r_destination_model_part = mrModel.GetModelPart(mDestinationModelPartName);

    KRATOS_ERROR_IF(r_source_model_part_nodes.size() !=
                    r_destination_model_part.Nodes().size())
        << "Number of nodes mismatch in sourc " << mSourceModelPartName
        << " and destination " << mDestinationModelPartName
        << " [ number of nodes in source = " << r_source_model_part_nodes.size()
        << ", number of nodes in destination model part = "
        << r_destination_model_part.Nodes().size() << " ].\n";

    block_for_each(r_source_model_part_nodes, [&](const NodeType& rSourceNode) {
        auto& r_destination_node = r_destination_model_part.GetNode(rSourceNode.Id());
        for (const auto& copy_variable_data : mCopyVariableDataList) {
            copy_variable_data.CopyData(rSourceNode, r_destination_node);
        }
    });

    // In here no mpi communication isrequired because we assume source model part is properly synchronized.

    if (mEchoLevel > 0) {
        std::stringstream msg;
        msg << "Following variable data was copied successfully from "
            << mSourceModelPartName << " to " << mDestinationModelPartName << ":\n";

        for (const auto& copy_variable_data : mCopyVariableDataList) {
            msg << "\t" << copy_variable_data.Info() << "\n";
        }

        KRATOS_INFO(this->Info()) << msg.str();
    }

    KRATOS_CATCH("");
}

void RansVariableDataTransferProcess::UpdateCopyExecutionPointsList(const std::vector<std::string>& rCopyExecutionPointsList)
{
    KRATOS_TRY

    mExecutionPointsList.clear();

    for (const auto& copy_execution_point : rCopyExecutionPointsList) {
        if (copy_execution_point == "initialize") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE);
        } else if (copy_execution_point == "initialize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE_SOLUTION_STEP);
        } else if (copy_execution_point == "before_coupling_solve_step") {
            mExecutionPointsList.push_back(ExecutionPoint::BEFORE_COUPLING_SOLVE_STEP);
        } else if (copy_execution_point == "execute") {
            mExecutionPointsList.push_back(ExecutionPoint::EXECUTE);
        } else if (copy_execution_point == "after_coupling_solve_step") {
            mExecutionPointsList.push_back(ExecutionPoint::AFTER_COUPLING_SOLVE_STEP);
        } else if (copy_execution_point == "finalize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE_SOLUTION_STEP);
        } else if (copy_execution_point == "finalize") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE);
        } else {
            KRATOS_ERROR << "Unsupported copy execution point provided. [ copy_execution_point = " << copy_execution_point << " ]. Supported points are: \n\n"
                        << "\tinitialize\n"
                        << "\tinitialize_solution_step\n"
                        << "\tbefore_coupling_solve_step\n"
                        << "\texecute\n"
                        << "\tafter_coupling_solve_step\n"
                        << "\tfinalize_solution_step\n"
                        << "\tfinalize\n";
        }
    }

    KRATOS_CATCH("");
}

void RansVariableDataTransferProcess::UpdateCopyVariableDataList(const std::vector<CopyVariableDataListItem>& rCopyVariableDataList)
{
    KRATOS_TRY

    mCopyVariableDataList.clear();

    for (const auto&  r_copy_variable_data_item : rCopyVariableDataList) {
        const auto& source_variable_name = std::get<0>(r_copy_variable_data_item);
        KRATOS_ERROR_IF_NOT(KratosComponents<VariableType>::Has(source_variable_name)) << "Source variable " << source_variable_name << " not found in double variables database." << std::endl;
        const VariableType& r_source_variable = KratosComponents<VariableType>::Get(source_variable_name);

        const bool is_source_variable_historical = std::get<1>(r_copy_variable_data_item);

        const auto& destination_variable_name = std::get<2>(r_copy_variable_data_item);
        KRATOS_ERROR_IF_NOT(KratosComponents<VariableType>::Has(destination_variable_name)) << "Destination variable " << destination_variable_name << " not found in double variables database." << std::endl;
        const VariableType& r_destination_variable = KratosComponents<VariableType>::Get(destination_variable_name);

        const bool is_destination_variable_historical = std::get<3>(r_copy_variable_data_item);

        mCopyVariableDataList.push_back(CopyVariableData(
            r_source_variable, is_source_variable_historical,
            r_destination_variable, is_destination_variable_historical));
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
