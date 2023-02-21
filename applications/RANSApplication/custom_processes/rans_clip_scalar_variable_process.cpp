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

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_variable_utilities.h"

// Include base h
#include "rans_clip_scalar_variable_process.h"

namespace Kratos
{
RansClipScalarVariableProcess::RansClipScalarVariableProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mVariableName = rParameters["variable_name"].GetString();
    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mMinValue = rParameters["min_value"].GetDouble();
    mMaxValue = rParameters["max_value"].GetDouble();

    KRATOS_CATCH("");
}

int RansClipScalarVariableProcess::Check()
{
    KRATOS_TRY

    const auto& r_scalar_variable = KratosComponents<Variable<double>>::Get(mVariableName);
    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_scalar_variable))
        << mVariableName << " is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansClipScalarVariableProcess::Execute()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_scalar_variable = KratosComponents<Variable<double>>::Get(mVariableName);

    unsigned int nodes_below, nodes_above;
    std::tie(nodes_below, nodes_above) = RansVariableUtilities::ClipScalarVariable(
        mMinValue, mMaxValue, r_scalar_variable, r_model_part);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && (nodes_below > 0 || nodes_above > 0))
        << mVariableName << " is clipped between [ " << mMinValue << ", "
        << mMaxValue << " ]. [ " << nodes_below << " nodes < " << mMinValue
        << " and " << nodes_above << " nodes > " << mMaxValue << " out of "
        << r_model_part.GetCommunicator().GlobalNumberOfNodes() << " total nodes in " << mModelPartName << " ].\n";

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

const Parameters RansClipScalarVariableProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "max_value"       : 1e+30
        })");

    return default_parameters;
}

} // namespace Kratos.
