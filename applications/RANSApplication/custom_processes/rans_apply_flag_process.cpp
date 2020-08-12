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
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "rans_apply_flag_process.h"

namespace Kratos
{
RansApplyFlagProcess::RansApplyFlagProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0,
            "flag_variable_name"             : "PLEASE_SPECIFY_FLAG_VARIABLE_NAME",
            "flag_variable_value"            : true,
            "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mFlagVariableName = mrParameters["flag_variable_name"].GetString();
    mFlagVariableValue = mrParameters["flag_variable_value"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartsForConditionFlags =
        mrParameters["apply_to_model_part_conditions"].GetStringArray();

    KRATOS_CATCH("");
}

void RansApplyFlagProcess::ExecuteInitialize()
{
    ApplyNodeFlags();

    if (mModelPartsForConditionFlags.size() == 1 &&
        mModelPartsForConditionFlags[0] == "ALL_MODEL_PARTS") {
        mModelPartsForConditionFlags.clear();
        for (const std::string& model_part_name : mrModel.GetModelPartNames()) {
            mModelPartsForConditionFlags.push_back(model_part_name);
        }
    }

    for (std::string model_part_name : mModelPartsForConditionFlags) {
        auto& r_model_part = mrModel.GetModelPart(model_part_name);
        ApplyConditionFlags(r_model_part);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << mFlagVariableName << " condition flag set for " << mModelPartName << ".\n";
}

std::string RansApplyFlagProcess::Info() const
{
    return std::string("RansApplyFlagProcess");
}

void RansApplyFlagProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansApplyFlagProcess::PrintData(std::ostream& rOStream) const
{
}

void RansApplyFlagProcess::ApplyNodeFlags()
{
    KRATOS_TRY

    auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();
    const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

    BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each(
        [&](ModelPart::NodeType& rNode) {
            rNode.Set(r_flag, mFlagVariableValue);
        });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << mFlagVariableName << " is set to nodes " << mFlagVariableValue
        << " in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

void RansApplyFlagProcess::ApplyConditionFlags(ModelPart& rModelPart)
{
    KRATOS_TRY

    auto& r_conditions = rModelPart.Conditions();
    const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

    BlockPartition<ModelPart::ConditionsContainerType>(r_conditions)
        .for_each([&](ModelPart::ConditionType& rCondition) {
            auto& r_condition_geometry = rCondition.GetGeometry();
            const int number_of_nodes = r_condition_geometry.PointsNumber();

            bool condition_flag = mFlagVariableValue;
            for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
                if (r_condition_geometry[i_node].Is(r_flag) != mFlagVariableValue) {
                    condition_flag = !condition_flag;
                    break;
                }
            }
            rCondition.Set(r_flag, condition_flag);
        });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << mFlagVariableName << " flags set for conditions in "
        << rModelPart.Name() << ".\n";

    KRATOS_CATCH("");
}

} // namespace Kratos.
