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
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_epsilon_turbulent_mixing_length_inlet_process.h"

namespace Kratos
{
RansEpsilonTurbulentMixingLengthInletProcess::RansEpsilonTurbulentMixingLengthInletProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mTurbulentMixingLength = rParameters["turbulent_mixing_length"].GetDouble();
    mIsConstrained = rParameters["is_fixed"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mTurbulentMixingLength < std::numeric_limits<double>::epsilon())
        << "turbulent_mixing_length should be greater than zero.\n";

    KRATOS_ERROR_IF(mMinValue < 0.0) << "Minimum turbulent energy dissipation "
                                        "rate needs to be positive in the "
                                        "modelpart "
                                     << mModelPartName << "\n.";

    KRATOS_CATCH("");
}

void RansEpsilonTurbulentMixingLengthInletProcess::ExecuteInitialize()
{
    if (mIsConstrained) {
        auto& r_model_part = mrModel.GetModelPart(mModelPartName);

        VariableUtils().ApplyFixity(TURBULENT_ENERGY_DISSIPATION_RATE, true,
                                    r_model_part.Nodes());

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed TURBULENT_ENERGY_DISSIPATION_RATE dofs in "
            << mModelPartName << ".\n";
    }
}

void RansEpsilonTurbulentMixingLengthInletProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double c_mu_75 = std::pow(r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU], 0.75);
    auto& r_nodes = r_model_part.Nodes();

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = std::max(
            c_mu_75 * std::pow(std::max(tke, 0.0), 1.5) / mTurbulentMixingLength, mMinValue);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied epsilon values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansEpsilonTurbulentMixingLengthInletProcess::Check()
{
    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;
}

std::string RansEpsilonTurbulentMixingLengthInletProcess::Info() const
{
    return std::string("RansEpsilonTurbulentMixingLengthInletProcess");
}

void RansEpsilonTurbulentMixingLengthInletProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansEpsilonTurbulentMixingLengthInletProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansEpsilonTurbulentMixingLengthInletProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_mixing_length" : 0.005,
            "echo_level"              : 0,
            "is_fixed"                : true,
            "min_value"               : 1e-14
        })");

    return default_parameters;
}

} // namespace Kratos.
