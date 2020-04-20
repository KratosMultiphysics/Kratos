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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_wall_function_update_process.h"

namespace Kratos
{
RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "von_karman"      : 0.41,
            "beta"            : 5.2
        })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mVonKarman = rParameters["von_karman"].GetDouble();
    mBeta = rParameters["beta"].GetDouble();

    KRATOS_CATCH("");
}

RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(Model& rModel,
                                                             const std::string& rModelPartName,
                                                             const double VonKarman,
                                                             const double Beta,
                                                             const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mVonKarman(VonKarman),
      mBeta(Beta),
      mEchoLevel(EchoLevel)
{
}

int RansWallFunctionUpdateProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VELOCITY);

    return 0;

    KRATOS_CATCH("");
}

void RansWallFunctionUpdateProcess::ExecuteInitialize()
{
    CalculateConditionNeighbourCount();
}

void RansWallFunctionUpdateProcess::CalculateConditionNeighbourCount()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        ConditionGeometryType& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour conditions in " << mModelPartName << ".\n";
}

void RansWallFunctionUpdateProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    ModelPart::NodesContainerType& r_nodes = r_model_part.Nodes();
    VariableUtils().SetHistoricalVariableToZero(RANS_Y_PLUS, r_nodes);
    VariableUtils().SetNonHistoricalVariableToZero(FRICTION_VELOCITY, r_nodes);

    ModelPart::ConditionsContainerType& r_conditions = r_model_part.Conditions();
    const int number_of_conditions = r_conditions.size();

#pragma omp parallel
    {
        array_1d<double, 3> r_friction_velocity;
#pragma omp for
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
        {
            ModelPart::ConditionType& r_condition = *(r_conditions.begin() + i_cond);

            const array_1d<double, 3>& r_wall_cell_center_velocity =
                RansCalculationUtilities::CalculateWallVelocity(r_condition);
            const double wall_cell_center_velocity_magnitude =
                norm_2(r_wall_cell_center_velocity);
            const array_1d<double, 3>& r_normal = r_condition.GetValue(NORMAL);

            const double wall_height =
                RansCalculationUtilities::CalculateWallHeight(r_condition, r_normal);

            const double nu = RansCalculationUtilities::EvaluateInParentCenter(
                KINEMATIC_VISCOSITY, r_condition);

            double y_plus{0.0}, u_tau{0.0};
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_cell_center_velocity_magnitude, wall_height,
                nu, mVonKarman, mBeta);

            if (wall_cell_center_velocity_magnitude > 0.0)
            {
                noalias(r_friction_velocity) =
                    r_wall_cell_center_velocity * (u_tau / wall_cell_center_velocity_magnitude);
            }
            else
            {
                noalias(r_friction_velocity) = ZeroVector(3);
            }

            r_condition.SetValue(RANS_Y_PLUS, y_plus);
            r_condition.SetValue(FRICTION_VELOCITY, r_friction_velocity);

            ModelPart::ConditionType::GeometryType& r_geometry =
                r_condition.GetGeometry();
            const int number_of_nodes = r_geometry.PointsNumber();
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                ModelPart::NodeType& r_node = r_geometry[i_node];
                const array_1d<double, 3>& r_u_tau = r_node.GetValue(FRICTION_VELOCITY);
                const double inv_number_of_neighbour_conditions =
                    1.0 / static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));

                r_node.SetLock();
                r_node.FastGetSolutionStepValue(RANS_Y_PLUS) +=
                    y_plus * inv_number_of_neighbour_conditions;
                r_node.SetValue(FRICTION_VELOCITY,
                                r_u_tau + r_friction_velocity * inv_number_of_neighbour_conditions);
                r_node.UnSetLock();
            }
        }
    }

    r_model_part.GetCommunicator().AssembleCurrentData(RANS_Y_PLUS);
    r_model_part.GetCommunicator().AssembleNonHistoricalData(FRICTION_VELOCITY);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated wall function based y_plus for " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansWallFunctionUpdateProcess::Info() const
{
    return std::string("RansWallFunctionUpdateProcess");
}

void RansWallFunctionUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansWallFunctionUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
