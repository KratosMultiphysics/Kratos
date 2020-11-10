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
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Include base h
#include "rans_apply_flag_to_skin_process.h"

namespace Kratos
{
RansApplyFlagToSkinProcess::RansApplyFlagToSkinProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mFlagVariableName = rParameters["flag_variable_name"].GetString();
    mFlagVariableValue = rParameters["flag_variable_value"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartsForConditionFlags = rParameters["apply_to_model_part_conditions"].GetStringArray();

    KRATOS_CATCH("");
}

void RansApplyFlagToSkinProcess::ExecuteInitialize()
{
    ApplyNodeFlags();

    if (mModelPartsForConditionFlags.size() == 1 &&
        mModelPartsForConditionFlags[0] == "ALL_MODEL_PARTS") {
        mModelPartsForConditionFlags.clear();
        for (const auto& model_part_name : mrModel.GetModelPartNames()) {
            mModelPartsForConditionFlags.push_back(model_part_name);
        }
    }

    for (const auto& model_part_name : mModelPartsForConditionFlags) {
        auto& r_model_part = mrModel.GetModelPart(model_part_name);
        ApplyConditionFlags(r_model_part);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << mFlagVariableName << " condition flag set for " << mModelPartName << ".\n";
}

std::string RansApplyFlagToSkinProcess::Info() const
{
    return std::string("RansApplyFlagToSkinProcess");
}

void RansApplyFlagToSkinProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansApplyFlagToSkinProcess::PrintData(std::ostream& rOStream) const
{
}

void RansApplyFlagToSkinProcess::ApplyNodeFlags()
{
    KRATOS_TRY

    auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();
    const auto& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

    VariableUtils().SetFlag(r_flag, mFlagVariableValue, r_nodes);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << mFlagVariableName << " is set to nodes " << mFlagVariableValue
        << " in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

void RansApplyFlagToSkinProcess::ApplyConditionFlags(
    ModelPart& rModelPart)
{
    KRATOS_TRY

    auto& r_conditions = rModelPart.Conditions();
    const auto& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

    block_for_each(r_conditions, [&](ModelPart::ConditionType& rCondition) {
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

const Parameters RansApplyFlagToSkinProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0,
            "flag_variable_name"             : "PLEASE_SPECIFY_FLAG_VARIABLE_NAME",
            "flag_variable_value"            : true,
            "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
        })");

    return default_parameters;
}

} // namespace Kratos.
