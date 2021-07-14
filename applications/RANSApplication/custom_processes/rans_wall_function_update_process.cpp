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
#include "custom_utilities/fluid_calculation_utilities.h"
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

    KRATOS_CATCH("");
}

RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mEchoLevel(EchoLevel)
{
}

int RansWallFunctionUpdateProcess::Check()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VELOCITY))
        << "VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    if (r_model_part.GetValue(RANS_IS_WALL_FUNCTION_ACTIVE) == 0) {
        return 0;
    }

    block_for_each(r_model_part.Conditions(), [&](const ModelPart::ConditionType& rCondition) {
        KRATOS_ERROR_IF_NOT(rCondition.Has(GAUSS_RANS_Y_PLUS))
            << "GAUSS_RANS_Y_PLUS is not found in condition data value container [ Condition.Id() = "
            << rCondition.Id() << " ].\n";

        const std::size_t number_of_gauss_points = rCondition.GetGeometry().IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_2);
        KRATOS_ERROR_IF(rCondition.GetValue(GAUSS_RANS_Y_PLUS).size() != number_of_gauss_points)
            << "GAUSS_RANS_Y_PLUS is not properly initialized. [ GAUSS_RANS_Y_PLUS.size() = " << rCondition.GetValue(GAUSS_RANS_Y_PLUS).size()
            << ", required size = " << number_of_gauss_points << " ].\n";
    });

    return 0;

    KRATOS_CATCH("");
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

void RansWallFunctionUpdateProcess::ExecuteAfterCouplingSolveStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (r_model_part.GetValue(RANS_IS_WALL_FUNCTION_ACTIVE) == 0) {
        return;
    }

    const auto& r_process_info = r_model_part.GetProcessInfo();
    const double von_karman = r_process_info[VON_KARMAN];

    auto& r_conditions = r_model_part.Conditions();

    block_for_each(r_conditions, std::pair<Vector, Matrix>(), [&](ModelPart::ConditionType& rCondition, std::pair<Vector, Matrix>& rTLS) {
        const auto& r_geometry = rCondition.GetGeometry();

        // get parent element
        auto& r_parent_element = r_geometry.GetValue(NEIGHBOUR_ELEMENTS)[0];

        // get fluid properties from parent element
        const auto& r_elem_properties = r_parent_element.GetProperties();
        const double rho = r_elem_properties[DENSITY];
        const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;

        // get surface properties
        const auto& r_cond_properties = rCondition.GetProperties();
        const double y_plus_limit = r_cond_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];
        const double beta = r_cond_properties[WALL_SMOOTHNESS_BETA];

        auto& gauss_weights = rTLS.first;
        auto& shape_functions = rTLS.second;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2,
            gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const double wall_height = rCondition.GetValue(DISTANCE);

        double tke;
        array_1d<double, 3> wall_velocity;

        auto& gauss_y_plus = rCondition.GetValue(GAUSS_RANS_Y_PLUS);

        for (size_t g = 0; g < num_gauss_points; ++g) {
            double& y_plus = gauss_y_plus[g];
            const auto& gauss_shape_functions = row(shape_functions, g);

            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, gauss_shape_functions,
                std::tie(wall_velocity, VELOCITY),
                std::tie(tke, TURBULENT_KINETIC_ENERGY));

            const double wall_velocity_magnitude = norm_2(wall_velocity);

            double u_tau{0.0};
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, wall_height, nu,
                von_karman, beta);
            y_plus = std::max(y_plus, y_plus_limit);
        }
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
            "echo_level"      : 0
        })");

    return default_parameters;
}

} // namespace Kratos.
