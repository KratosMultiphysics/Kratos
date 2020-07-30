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

#include "custom_utilities/rans_check_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

// Include base h
#include "rans_clip_scalar_variable_process.h"

namespace Kratos
{
RansClipScalarVariableProcess::RansClipScalarVariableProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "max_value"       : 1e+30
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = mrParameters["variable_name"].GetString();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mMinValue = mrParameters["min_value"].GetDouble();
    mMaxValue = mrParameters["max_value"].GetDouble();

    KRATOS_CATCH("");
}

int RansClipScalarVariableProcess::Check()
{
    KRATOS_TRY

    const Variable<double>& r_scalar_variable =
        KratosComponents<Variable<double>>::Get(mVariableName);


    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        mrModel.GetModelPart(mModelPartName), r_scalar_variable);

    return 0;

    KRATOS_CATCH("");
}

void RansClipScalarVariableProcess::Execute()
{
    KRATOS_TRY

    const Variable<double>& r_scalar_variable =
        KratosComponents<Variable<double>>::Get(mVariableName);

    unsigned int nodes_below, nodes_above, total_nodes;

    RansVariableUtilities::ClipScalarVariable(
        nodes_below, nodes_above, total_nodes, mMinValue, mMaxValue,
        r_scalar_variable, mrModel.GetModelPart(mModelPartName));

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && (nodes_below > 0 || nodes_above > 0))
        << mVariableName << " is clipped between [ " << mMinValue << ", "
        << mMaxValue << " ]. [ " << nodes_below << " nodes < " << mMinValue
        << " and " << nodes_above << " nodes > " << mMaxValue << " out of "
        << total_nodes << " total nodes in " << mModelPartName << " ].\n";

    KRATOS_CATCH("");
}

std::string RansClipScalarVariableProcess::Info() const
{
    return std::string("RansClipScalarVariableProcess");
}

void RansClipScalarVariableProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansClipScalarVariableProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
