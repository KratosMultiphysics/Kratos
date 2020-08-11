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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_wall_function_update_process.h"

namespace Kratos
{
RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "c_mu"            : 0.09
        })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mVonKarman = rParameters["von_karman"].GetDouble();
    mBeta = rParameters["beta"].GetDouble();
    mCmu = rParameters["c_mu"].GetDouble();

    KRATOS_CATCH("");
}

RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(
    Model& rModel,
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

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(KINEMATIC_VISCOSITY))
        << "KINEMATIC_VISCOSITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VELOCITY))
        << "VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansWallFunctionUpdateProcess::ExecuteInitialize()
{
    CalculateConditionNeighbourCount();
}

void RansWallFunctionUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        this->Execute();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansWallFunctionUpdateProcess::CalculateConditionNeighbourCount()
{
    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
        auto& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        auto& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            auto& r_node = r_geometry[i_node];
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

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double c_mu_25 = std::pow(mCmu, 0.25);
    const double y_plus_limit =
        r_model_part.GetProcessInfo()[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    auto& r_conditions = r_model_part.Conditions();
    const int number_of_conditions = r_conditions.size();

#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
        auto& r_condition = *(r_conditions.begin() + i_cond);

        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_condition.GetGeometry(), r_condition.GetIntegrationMethod(),
            gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const array_1d<double, 3>& r_normal = r_condition.GetValue(NORMAL);
        const double wall_height =
            RansCalculationUtilities::CalculateWallHeight(r_condition, r_normal);

        double condition_y_plus{0.0};
        array_1d<double, 3> condition_u_tau = ZeroVector(3);

        for (size_t g = 0; g < num_gauss_points; ++g) {
            const Vector& gauss_shape_functions = row(shape_functions, g);
            const array_1d<double, 3>& r_wall_velocity =
                RansCalculationUtilities::EvaluateInPoint(
                    r_condition.GetGeometry(), VELOCITY, gauss_shape_functions);
            const double wall_velocity_magnitude = norm_2(r_wall_velocity);
            const double nu = RansCalculationUtilities::EvaluateInPoint(
                r_condition.GetGeometry(), KINEMATIC_VISCOSITY, gauss_shape_functions);
            const double tke = RansCalculationUtilities::EvaluateInPoint(
                r_condition.GetGeometry(), TURBULENT_KINETIC_ENERGY, gauss_shape_functions);

            double y_plus{0.0}, u_tau{0.0};
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, wall_height, nu, mVonKarman, mBeta);
            y_plus = std::max(y_plus, y_plus_limit);

            u_tau = RansCalculationUtilities::SoftMax(
                c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                wall_velocity_magnitude / (std::log(y_plus) / mVonKarman + mBeta));

            condition_y_plus += y_plus;

            if (wall_velocity_magnitude > 0.0) {
                noalias(condition_u_tau) += r_wall_velocity * u_tau / wall_velocity_magnitude;
            }
        }

        const double inv_number_of_gauss_points =
            static_cast<double>(1.0 / num_gauss_points);

        r_condition.SetValue(RANS_Y_PLUS, condition_y_plus * inv_number_of_gauss_points);
        r_condition.SetValue(FRICTION_VELOCITY, condition_u_tau * inv_number_of_gauss_points);
    }

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
