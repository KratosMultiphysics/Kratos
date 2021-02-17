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
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_k_turbulent_intensity_inlet_process.h"

namespace Kratos
{
RansKTurbulentIntensityInletProcess::RansKTurbulentIntensityInletProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mTurbulentIntensity = rParameters["turbulent_intensity"].GetDouble();
    mIsConstrained = rParameters["is_fixed"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mTurbulentIntensity < 0.0)
        << "Turbulent intensity needs to be positive in the modelpart "
        << mModelPartName << "\n.";
    KRATOS_ERROR_IF(mMinValue < 0.0)
        << "Minimum turbulent kinetic energy needs to be positive in the "
           "modelpart "
        << mModelPartName << "\n.";

    KRATOS_CATCH("");
}

void RansKTurbulentIntensityInletProcess::ExecuteInitialize()
{
    if (mIsConstrained) {
        auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();

        VariableUtils().ApplyFixity(TURBULENT_KINETIC_ENERGY, true,
                                    r_nodes);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed TURBULENT_KINETIC_ENERGY dofs in " << mModelPartName << ".\n";
    }
}

void RansKTurbulentIntensityInletProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        const auto& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        const double velocity_magnitude = norm_2(r_velocity);

        rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = std::max(
            1.5 * std::pow(mTurbulentIntensity * velocity_magnitude, 2), mMinValue);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied k values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansKTurbulentIntensityInletProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VELOCITY))
        << "VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

std::string RansKTurbulentIntensityInletProcess::Info() const
{
    return std::string("RansKTurbulentIntensityInletProcess");
}

void RansKTurbulentIntensityInletProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansKTurbulentIntensityInletProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansKTurbulentIntensityInletProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_intensity" : 0.05,
            "echo_level"          : 0,
            "is_fixed"            : true,
            "min_value"           : 1e-14
        })");

    return default_parameters;
}

} // namespace Kratos.
