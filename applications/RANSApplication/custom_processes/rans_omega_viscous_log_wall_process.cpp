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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_omega_viscous_log_wall_process.h"
namespace Kratos
{
/// Constructor
RansOmegaViscousLogWallProcess::RansOmegaViscousLogWallProcess(
    Model& rModel,
    Parameters& rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mIsConstrained = rParameters["is_fixed"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mCalculationStepIndex = rParameters["calculation_step_index"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mMinValue < 0.0) << "Minimum turbulent energy dissipation "
                                        "rate needs to be positive in the "
                                        "modelpart "
                                     << mModelPartName << "\n.";

    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

void RansOmegaViscousLogWallProcess::ExecuteInitialize()
{
    if (mIsConstrained) {
        auto& r_nodes = mrModel.GetModelPart(mModelPartName).Nodes();

        VariableUtils().ApplyFixity(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, true, r_nodes);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Fixed TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE dofs in "
            << mModelPartName << ".\n";
    }

    BaseType::ExecuteInitialize();
}

void RansOmegaViscousLogWallProcess::ExecuteOperation()
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Matrix, Geometry<Node>::ShapeFunctionsGradientsType>;

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const auto& r_process_info = r_model_part.GetProcessInfo();
    const double kappa = r_process_info[VON_KARMAN];
    const double c_mu_25 = std::pow(r_process_info[TURBULENCE_RANS_C_MU], 0.25);
    auto& r_conditions = r_model_part.Conditions();

    VariableUtils().SetHistoricalVariableToZero(
        TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_model_part.Nodes());

    block_for_each(r_conditions, tls_type(), [&](ModelPart::ConditionType& rCondition, tls_type& rTLS) {
        auto& Ws = std::get<0>(rTLS);
        auto& Ns = std::get<1>(rTLS);
        auto& dNdXs = std::get<2>(rTLS);

        // get parent element
        const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
        const auto& r_parent_element_geometry = r_parent_element.GetGeometry();

        // get fluid properties from parent element
        const auto& r_elem_properties = r_parent_element.GetProperties();
        const double rho = r_elem_properties[DENSITY];
        const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;
        const double y = rCondition.GetValue(DISTANCE);
        array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
        normal /= norm_2(normal);

        // get surface properties from condition
        const auto& r_cond_properties = rCondition.GetProperties();
        const double beta = r_cond_properties.GetValue(WALL_SMOOTHNESS_BETA);
        const double y_plus_limit = r_cond_properties.GetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT);
        const double k_s = r_cond_properties.GetValue(WALL_ROUGHNESS_HEIGHT);

        // Get Shape function data
        RansCalculationUtilities::CalculateGeometryData(
            r_parent_element_geometry,
            GeometryData::IntegrationMethod::GI_GAUSS_1, Ws, Ns, dNdXs);

        array_1d<double, 3> wall_velocity, fluid_velocity, mesh_velocity;
        double tke;
        FluidCalculationUtilities::EvaluateInPoint(
            r_parent_element_geometry, row(Ns, 0), mCalculationStepIndex,
            std::tie(tke, TURBULENT_KINETIC_ENERGY),
            std::tie(fluid_velocity, VELOCITY),
            std::tie(mesh_velocity, MESH_VELOCITY));

        noalias(wall_velocity) = fluid_velocity - mesh_velocity;
        const double wall_velocity_magnitude =
            std::sqrt(std::pow(norm_2(wall_velocity), 2) -
                      std::pow(inner_prod(wall_velocity, normal), 2));

        double u_tau, y_plus;
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, y, nu, kappa, beta);

        double omega;
        if (y_plus > y_plus_limit) {
            // log region
            const double tke_u_tau = c_mu_25 * std::sqrt(std::max(tke, 0.0));
            omega = std::max(u_tau, tke_u_tau) / (c_mu_25 * kappa * y);
        } else {
            // linear region

            const double k_s_plus = std::max(1.0, k_s * u_tau / nu);

            double omega_wall_plus = 100.0 / k_s_plus;
            if (k_s_plus < 25.0) {
                omega_wall_plus = std::pow(50.0 / k_s_plus, 2);
            }

            omega = (std::pow(u_tau, 2) / nu) *
                    std::min(omega_wall_plus, 6.0 / (0.075 * std::pow(y_plus, 2)));
        }

        auto& r_condition_geometry = rCondition.GetGeometry();
        omega = std::max(mMinValue, omega / r_condition_geometry.PointsNumber());
        for (auto& r_node : r_condition_geometry) {
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) += omega;
            r_node.UnSetLock();
        }
    });

    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied omega values to " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansOmegaViscousLogWallProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VELOCITY))
        << "VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(MESH_VELOCITY))
        << "MESH_VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansOmegaViscousLogWallProcess::Info() const
{
    return std::string("RansOmegaViscousLogWallProcess");
}

/// Print information about this object.
void RansOmegaViscousLogWallProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansOmegaViscousLogWallProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansOmegaViscousLogWallProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"              : 0,
        "is_fixed"                : true,
        "min_value"               : 1e-12,
        "execution_points"        : ["initialize_solution_step"],
        "calculation_step_index"  : 1
    })");

    return default_parameters;
}

} // namespace Kratos.