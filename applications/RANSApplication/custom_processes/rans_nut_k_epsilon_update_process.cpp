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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_k_epsilon_update_process.h"

namespace Kratos
{
RansNutKEpsilonUpdateProcess::RansNutKEpsilonUpdateProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mCmu = rParameters["c_mu"].GetDouble();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_CATCH("");
}

RansNutKEpsilonUpdateProcess::RansNutKEpsilonUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double Cmu,
    const double MinValue,
    const int EchoLevel)
: mrModel(rModel),
  mModelPartName(rModelPartName),
  mCmu(Cmu),
  mMinValue(MinValue),
  mEchoLevel(EchoLevel)
{
}

int RansNutKEpsilonUpdateProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_VISCOSITY))
        << "TURBULENT_VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansNutKEpsilonUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        this->ExecuteAfterCouplingSolveStep();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansNutKEpsilonUpdateProcess::ExecuteAfterCouplingSolveStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_nodes = r_model_part.Nodes();

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        const double epsilon =
            rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

        double& nu_t = rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY);

        if (epsilon > 0.0) {
            const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            nu_t = mCmu * std::pow(tke, 2) / epsilon;
        } else {
            nu_t = mMinValue;
        }

        rNode.FastGetSolutionStepValue(VISCOSITY) =
            rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) + nu_t;
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t for nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansNutKEpsilonUpdateProcess::Info() const
{
    return std::string("RansNutKEpsilonUpdateProcess");
}

void RansNutKEpsilonUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutKEpsilonUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansNutKEpsilonUpdateProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "min_value"       : 1e-15
        })");
    return default_parameters;
}

} // namespace Kratos.
