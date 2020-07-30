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
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

#include "custom_utilities/rans_check_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

// Include base h
#include "rans_check_scalar_bounds_process.h"

namespace Kratos
{
RansCheckScalarBoundsProcess::RansCheckScalarBoundsProcess(Model& rModel, Parameters rParameters)
    : Process(), mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE"
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = mrParameters["variable_name"].GetString();
    mModelPartName = mrParameters["model_part_name"].GetString();

    KRATOS_CATCH("");
}

int RansCheckScalarBoundsProcess::Check()
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

void RansCheckScalarBoundsProcess::Execute()
{
    KRATOS_TRY

    const Variable<double>& r_scalar_variable =
        KratosComponents<Variable<double>>::Get(mVariableName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double min_value =
        RansVariableUtilities::GetMinimumScalarValue(r_model_part, r_scalar_variable);
    const double max_value =
        RansVariableUtilities::GetMaximumScalarValue(r_model_part, r_scalar_variable);

    KRATOS_INFO(this->Info())
        << r_scalar_variable.Name() << " is bounded between [ " << min_value
        << ", " << max_value << " ] in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansCheckScalarBoundsProcess::Info() const
{
    return std::string("RansCheckScalarBoundsProcess");
}

void RansCheckScalarBoundsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansCheckScalarBoundsProcess::PrintData(std::ostream& rOStream) const
{
}
} // namespace Kratos.
