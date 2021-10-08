//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//                   Bence Rochlitz
//
//  Supervised by:   Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_omega_automatic_inlet_process.h"
namespace Kratos
{
/// Constructor
RansOmegaAutomaticInletProcess::RansOmegaAutomaticInletProcess(
    Model& rModel,
    Parameters& rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mTurbulentMixingLength = rParameters["turbulent_mixing_length"].GetDouble();
    mKinematicViscosity = rParameters["kinematic_viscosity"].GetDouble();
    mBeta = rParameters["beta"].GetDouble();
    mWallLocation = rParameters["wall_location"].GetVector();
    mWallOutwardPointintUnitNormal = rParameters["wall_normal"].GetVector();
    mWallOutwardPointintUnitNormal = mWallOutwardPointintUnitNormal / norm_2(mWallOutwardPointintUnitNormal);
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

void RansOmegaAutomaticInletProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    const double c_mu_25 = std::pow(r_process_info[TURBULENCE_RANS_C_MU], 0.25);
    const double kappa = r_process_info[VON_KARMAN];
    auto& r_nodes = r_model_part.Nodes();

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        const double y = inner_prod(mWallLocation - rNode.Coordinates(), mWallOutwardPointintUnitNormal);
        if (y > 1e-12) {
            const double tke_sqrt = std::sqrt(std::max(rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY), 0.0));
            const double omega_turbulent_mixing_length =  tke_sqrt / (c_mu_25 * mTurbulentMixingLength);
            const double omega_vis = 6.0 * mKinematicViscosity / (mBeta * y * y);
            const double omega_log = tke_sqrt / (c_mu_25 * kappa * y);

            double& omega = rNode.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
            omega = std::max(std::sqrt(std::pow(omega_vis, 2) + std::pow(omega_log, 2) +
                                       std::pow(omega_turbulent_mixing_length, 2)),
                             mMinValue);
            rNode.Fix(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        } else {
            rNode.Free(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        }
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied omega values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansOmegaAutomaticInletProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansOmegaAutomaticInletProcess::Info() const
{
    return std::string("RansOmegaAutomaticInletProcess");
}

/// Print information about this object.
void RansOmegaAutomaticInletProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansOmegaAutomaticInletProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansOmegaAutomaticInletProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "turbulent_mixing_length" : 0.005,
        "echo_level"              : 0,
        "beta"                    : 0.075,
        "kinematic_viscosity"     : 1e-5,
        "wall_location"           : [0.0, 0.0, 0.0],
        "wall_normal"             : [1.0, 0.0, 0.0],
        "min_value"               : 1e-12
    })");

    return default_parameters;
}

} // namespace Kratos.