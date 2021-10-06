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
#include "rans_omega_viscous_log_binomial_wall_process.h"
namespace Kratos
{
/// Constructor
RansOmegaViscousLogBinomialWallProcess::RansOmegaViscousLogBinomialWallProcess(
    Model& rModel,
    Parameters& rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mIsConstrained = rParameters["is_fixed"].GetBool();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_ERROR_IF(mMinValue < 0.0) << "Minimum turbulent energy dissipation "
                                        "rate needs to be positive in the "
                                        "modelpart "
                                     << mModelPartName << "\n.";

    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

void RansOmegaViscousLogBinomialWallProcess::ExecuteInitialize()
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

void RansOmegaViscousLogBinomialWallProcess::ExecuteOperation()
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Matrix, Geometry<Node<3>>::ShapeFunctionsGradientsType>;

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const auto& r_process_info = r_model_part.GetProcessInfo();
    const double kappa = r_process_info[VON_KARMAN];
    const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];
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
        const array_1d<double, 3>& normal = rCondition.GetValue(NORMAL);

        // Get Shape function data
        RansCalculationUtilities::CalculateGeometryData(
            r_parent_element_geometry,
            GeometryData::IntegrationMethod::GI_GAUSS_1, Ws, Ns, dNdXs);

        const double y = RansCalculationUtilities::CalculateWallHeight(rCondition, normal);

        double tke;
        FluidCalculationUtilities::EvaluateInPoint(
            r_parent_element_geometry, row(Ns, 0),
            std::tie(tke, TURBULENT_KINETIC_ENERGY));

        const double omega_vis = 6.0 * nu / (0.075 * y * y);
        const double omega_log = std::sqrt(std::max(tke, 0.0)) / (c_mu * kappa * y);
        double omega = std::sqrt(omega_vis * omega_vis + omega_log * omega_log);

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

int RansOmegaViscousLogBinomialWallProcess::Check()
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
std::string RansOmegaViscousLogBinomialWallProcess::Info() const
{
    return std::string("RansOmegaViscousLogBinomialWallProcess");
}

/// Print information about this object.
void RansOmegaViscousLogBinomialWallProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansOmegaViscousLogBinomialWallProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansOmegaViscousLogBinomialWallProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"              : 0,
        "is_fixed"                : true,
        "min_value"               : 1e-12,
        "execution_points"        : ["initialize_solution_step"]
    })");

    return default_parameters;
}

} // namespace Kratos.