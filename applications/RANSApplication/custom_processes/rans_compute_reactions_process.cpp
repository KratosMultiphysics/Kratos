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
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_compute_reactions_process.h"

namespace Kratos
{
RansComputeReactionsProcess::RansComputeReactionsProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mPeriodic = rParameters["consider_periodic"].GetBool();

    KRATOS_CATCH("");
}

void RansComputeReactionsProcess::ExecuteInitialize()
{
    KRATOS_TRY

    if (mPeriodic)
    {
        this->CorrectPeriodicNodes(mrModel.GetModelPart(mModelPartName), NORMAL);
    }

    KRATOS_CATCH("");
}

void RansComputeReactionsProcess::CorrectPeriodicNodes(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable)
{
    KRATOS_TRY

    auto& r_nodes = rModelPart.Nodes();

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        if (rNode.Is(PERIODIC)) {
            const int slave_id = rNode.FastGetSolutionStepValue(PATCH_INDEX);
            const int master_id = rNode.Id();
            if (master_id < slave_id) {
                auto& master_value = rNode.FastGetSolutionStepValue(rVariable);

                auto& slave_node = rModelPart.GetNode(slave_id);
                auto& slave_value = slave_node.FastGetSolutionStepValue(rVariable);

                const double master_value_norm = norm_2(master_value);
                const double slave_value_norm = norm_2(slave_value);
                const double value_norm = master_value_norm + slave_value_norm;

                if (master_value_norm > 0.0) {
                    noalias(master_value) = master_value * (value_norm / master_value_norm);
                }

                if (slave_value_norm > 0.0) {
                    noalias(slave_value) = slave_value * (value_norm / slave_value_norm);
                }
            }
        }
    });

    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);

    KRATOS_CATCH("");
}

void RansComputeReactionsProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_nodes = r_model_part.Nodes();
    VariableUtils().SetHistoricalVariableToZero(REACTION, r_nodes);

    auto& r_conditions = r_model_part.Conditions();

    block_for_each(r_conditions, [&](ModelPart::ConditionType& rCondition) {
            CalculateReactionValues(rCondition);
        });

    r_model_part.GetCommunicator().AssembleCurrentData(REACTION);
    this->CorrectPeriodicNodes(r_model_part, REACTION);

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        const auto& pressure_force = rNode.FastGetSolutionStepValue(NORMAL) *
                                     (-1.0 * rNode.FastGetSolutionStepValue(PRESSURE));
        rNode.FastGetSolutionStepValue(REACTION) += pressure_force;
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Computed Reaction terms for " << mModelPartName << " for slip condition.\n";

    KRATOS_CATCH("");
}

int RansComputeReactionsProcess::Check()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(REACTION))
        << "REACTION is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(DENSITY))
        << "DENSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(PRESSURE))
        << "PRESSURE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(NORMAL))
        << "NORMAL is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

std::string RansComputeReactionsProcess::Info() const
{
    return std::string("RansComputeReactionsProcess");
}

void RansComputeReactionsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansComputeReactionsProcess::PrintData(std::ostream& rOStream) const
{
}

void RansComputeReactionsProcess::CalculateReactionValues(
    ModelPart::ConditionType& rCondition)
{
    const auto& r_friction_velocity = rCondition.GetValue(FRICTION_VELOCITY);
    const double u_tau = norm_2(r_friction_velocity);
    auto& r_geometry = rCondition.GetGeometry();

    const IndexType number_of_nodes = r_geometry.PointsNumber();

    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        auto& r_node = r_geometry[i_node];
        const double rho = r_node.FastGetSolutionStepValue(DENSITY);

        if (u_tau > 0.0)
        {
            const double shear_force = rho * std::pow(u_tau, 2) *
                                       r_geometry.DomainSize() /
                                       static_cast<double>(number_of_nodes);
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(REACTION) +=
                r_friction_velocity * (shear_force / u_tau);
            r_node.UnSetLock();
        }
    }
}

const Parameters RansComputeReactionsProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"              : 0,
            "consider_periodic"       : false
        })");

    return default_parameters;
}

} // namespace Kratos.