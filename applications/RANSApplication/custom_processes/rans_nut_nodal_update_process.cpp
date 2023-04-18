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
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
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

    return 0;

    KRATOS_CATCH("");
}

void RansNutNodalUpdateProcess::ExecuteInitialize()
{
    RansCalculationUtilities::CalculateNumberOfNeighbourEntities<ModelPart::ElementsContainerType>(
        mrModel.GetModelPart(mModelPartName), NUMBER_OF_NEIGHBOUR_ELEMENTS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour elements in " << mModelPartName << ".\n";
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

    // clear all the nodes
    VariableUtils().SetHistoricalVariableToZero(VISCOSITY, r_nodes);

    // compute nu with nu_t for elemental nodes
    block_for_each(r_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
        const double nu_t = rElement.GetValue(TURBULENT_VISCOSITY);
        for (auto& r_node : rElement.GetGeometry()) {
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(VISCOSITY) += (nu + nu_t) / r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            r_node.UnSetLock();
        }
    });

    // assemble elemental nu_t computations
    r_model_part.GetCommunicator().AssembleCurrentData(VISCOSITY);

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
