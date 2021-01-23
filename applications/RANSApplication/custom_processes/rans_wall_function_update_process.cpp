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
#include "utilities/parallel_utilities.h"

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

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

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
        this->ExecuteAfterCouplingSolveStep();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansWallFunctionUpdateProcess::CalculateConditionNeighbourCount()
{
    RansCalculationUtilities::CalculateNumberOfNeighbourEntities<ModelPart::ConditionsContainerType>(
        mrModel.GetModelPart(mModelPartName), NUMBER_OF_NEIGHBOUR_CONDITIONS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour conditions in " << mModelPartName << ".\n";
}

void RansWallFunctionUpdateProcess::ExecuteAfterCouplingSolveStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double c_mu_25 = std::pow(mCmu, 0.25);
    const double y_plus_limit =
        r_model_part.GetProcessInfo()[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    auto& r_conditions = r_model_part.Conditions();

    // TODO: Make vectors and matrices TLS
    block_for_each(r_conditions, [&](ModelPart::ConditionType& rCondition) {
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            rCondition.GetGeometry(), rCondition.GetIntegrationMethod(),
            gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const auto& r_normal = rCondition.GetValue(NORMAL);
        const double wall_height =
            RansCalculationUtilities::CalculateWallHeight(rCondition, r_normal);

        double condition_y_plus{0.0}, nu, tke;
        array_1d<double, 3> condition_u_tau = ZeroVector(3);
        array_1d<double, 3> wall_velocity;

        for (size_t g = 0; g < num_gauss_points; ++g) {
            const auto& gauss_shape_functions = row(shape_functions, g);

            RansCalculationUtilities::EvaluateInPoint(
                rCondition.GetGeometry(), gauss_shape_functions,
                std::tie(wall_velocity, VELOCITY), std::tie(nu, KINEMATIC_VISCOSITY),
                std::tie(tke, TURBULENT_KINETIC_ENERGY));

            const double wall_velocity_magnitude = norm_2(wall_velocity);

            double y_plus{0.0}, u_tau{0.0};
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, wall_height, nu, mVonKarman, mBeta);
            y_plus = std::max(y_plus, y_plus_limit);

            u_tau = RansCalculationUtilities::SoftMax(
                c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                wall_velocity_magnitude / (std::log(y_plus) / mVonKarman + mBeta));

            condition_y_plus += y_plus;

            if (wall_velocity_magnitude > 0.0) {
                noalias(condition_u_tau) += wall_velocity * u_tau / wall_velocity_magnitude;
            }
        }

        const double inv_number_of_gauss_points = 1.0 / num_gauss_points;

        rCondition.SetValue(RANS_Y_PLUS, condition_y_plus * inv_number_of_gauss_points);
        rCondition.SetValue(FRICTION_VELOCITY, condition_u_tau * inv_number_of_gauss_points);
    });

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

const Parameters RansWallFunctionUpdateProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "c_mu"            : 0.09
        })");

    return default_parameters;
}

} // namespace Kratos.
