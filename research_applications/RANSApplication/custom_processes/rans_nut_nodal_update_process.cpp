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
#include "rans_nut_nodal_update_process.h"

namespace Kratos
{
RansNutNodalUpdateProcess::RansNutNodalUpdateProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();

    KRATOS_CATCH("");
}

RansNutNodalUpdateProcess::RansNutNodalUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const int EchoLevel)
: mrModel(rModel),
  mModelPartName(rModelPartName),
  mEchoLevel(EchoLevel)
{
}

int RansNutNodalUpdateProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VISCOSITY))
        << "VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_VISCOSITY))
        << "TURBULENT_VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansNutNodalUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        this->ExecuteAfterCouplingSolveStep();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansNutNodalUpdateProcess::ExecuteAfterCouplingSolveStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_nodes = r_model_part.Nodes();

    const auto& r_properties = r_model_part.Elements().front().GetProperties();
    const double rho = r_properties.GetValue(DENSITY);
    const double nu = r_properties.GetValue(DYNAMIC_VISCOSITY) / rho;

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        rNode.FastGetSolutionStepValue(VISCOSITY) =
            nu + rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Updated nu_t for nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansNutNodalUpdateProcess::Info() const
{
    return std::string("RansNutNodalUpdateProcess");
}

void RansNutNodalUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutNodalUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansNutNodalUpdateProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level": 0
        })");
    return default_parameters;
}

} // namespace Kratos.
