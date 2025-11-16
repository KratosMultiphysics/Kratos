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

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_compute_y_plus_process.h"
namespace Kratos
{
/// Constructor
RansComputeYPlusProcess::RansComputeYPlusProcess(
    Model& rModel,
    Parameters& rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mOutputVariableName = rParameters["output_variable_name"].GetString();
    mKappa = rParameters["kappa"].GetDouble();
    mBeta = rParameters["beta"].GetDouble();

    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

void RansComputeYPlusProcess::ExecuteOperation()
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Matrix, Geometry<Node>::ShapeFunctionsGradientsType>;

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_output_variable = KratosComponents<Variable<double>>::Get(mOutputVariableName);

    block_for_each(r_model_part.Conditions(), tls_type(), [&](ModelPart::ConditionType& rCondition, tls_type& rTLS) {
        auto& Ws = std::get<0>(rTLS);
        auto& Ns = std::get<1>(rTLS);
        auto& dNdXs = std::get<2>(rTLS);

        // get parent element
        auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
        auto& r_parent_element_geometry = r_parent_element.GetGeometry();

        // get fluid properties from parent element
        const auto& r_elem_properties = r_parent_element.GetProperties();
        const double rho = r_elem_properties[DENSITY];
        const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;
        const double y = rCondition.GetValue(DISTANCE);
        array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
        normal /= norm_2(normal);

        // Get Shape function data
        RansCalculationUtilities::CalculateGeometryData(
            r_parent_element_geometry,
            GeometryData::IntegrationMethod::GI_GAUSS_1, Ws, Ns, dNdXs);

        array_1d<double, 3> wall_velocity, fluid_velocity, mesh_velocity;
        FluidCalculationUtilities::EvaluateInPoint(
            r_parent_element_geometry, row(Ns, 0),
            std::tie(fluid_velocity, VELOCITY),
            std::tie(mesh_velocity, MESH_VELOCITY));

        noalias(wall_velocity) = fluid_velocity - mesh_velocity;
        const double wall_velocity_magnitude =
            std::sqrt(std::pow(norm_2(wall_velocity), 2) -
                      std::pow(inner_prod(wall_velocity, normal), 2));

        double u_tau, y_plus;
        RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, y, nu, mKappa, mBeta);

        r_parent_element.SetValue(r_output_variable, y_plus);

    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Computed y+ values in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

int RansComputeYPlusProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(VELOCITY))
        << "VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(MESH_VELOCITY))
        << "MESH_VELOCITY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansComputeYPlusProcess::Info() const
{
    return std::string("RansComputeYPlusProcess");
}

/// Print information about this object.
void RansComputeYPlusProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansComputeYPlusProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansComputeYPlusProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
    {
        "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"              : 0,
        "kappa"                   : 0.41,
        "beta"                    : 5.2,
        "output_variable_name"    : "RANS_Y_PLUS",
        "execution_points"        : ["initialize_solution_step"]
    })");

    return default_parameters;
}

} // namespace Kratos.