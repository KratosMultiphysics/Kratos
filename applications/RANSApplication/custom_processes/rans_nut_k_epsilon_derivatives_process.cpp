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
#include <tuple>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/node.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_nut_k_epsilon_derivatives_process.h"

namespace Kratos
{
RansNutKEpsilonDerivativesProcess::RansNutKEpsilonDerivativesProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();

    KRATOS_CATCH("");
}

RansNutKEpsilonDerivativesProcess::RansNutKEpsilonDerivativesProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double MinValue,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mMinValue(MinValue),
      mEchoLevel(EchoLevel)
{
}

int RansNutKEpsilonDerivativesProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY))
        << "TURBULENT_KINETIC_ENERGY is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE))
        << "TURBULENT_ENERGY_DISSIPATION_RATE is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansNutKEpsilonDerivativesProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using ElementType = ModelPart::ElementType;

    using tls_type = std::tuple<Vector, Matrix, ShapeFunctionDerivativesArrayType, Vector>;

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_properties = r_model_part.GetProcessInfo();
    const double c_mu = r_properties[TURBULENCE_RANS_C_MU];
    const int derivatives_size = 2;

    // setting element nu_t values
    block_for_each(r_model_part.Elements(), tls_type(), [&](ElementType& rElement, tls_type& rTLS) {
        auto& Ws = std::get<0>(rTLS);
        auto& Ns = std::get<1>(rTLS);
        auto& dNdXs = std::get<2>(rTLS);
        auto& nu_t_derivatives = std::get<3>(rTLS);

        // reserve space for gauss nu_t derivatives w.r.t. gauss variables
        // VELOCITY_X, VELOCITY_Y, (VELOCITY_Z), TURBULENT_KINETIC_ENERGY, TURBULENT_ENERGY_DISSIPATION_RATE
        if (nu_t_derivatives.size() != derivatives_size) {
            nu_t_derivatives.resize(derivatives_size);
        }

        noalias(nu_t_derivatives) = ZeroVector(derivatives_size);

        // here the gauss integration method is fixed because,
        // the primal elements can have different gauss point integration methods
        // thus making it difficult to decide on which GP integration to use
        const auto& r_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;

        RansCalculationUtilities::CalculateGeometryData(
            rElement.GetGeometry(), r_integration_method, Ws, Ns, dNdXs);

        const Vector& N = row(Ns, 0);

        double tke, epsilon;
        FluidCalculationUtilities::EvaluateInPoint(rElement.GetGeometry(), N,
            std::tie(tke, TURBULENT_KINETIC_ENERGY),
            std::tie(epsilon, TURBULENT_ENERGY_DISSIPATION_RATE));

        double nu_t = 0.0;
        if (epsilon > 0.0) {
            nu_t = c_mu * std::pow(tke, 2) / epsilon;

            // compute derivative w.r.t. TURBULENT_KINETIC_ENERGY
            nu_t_derivatives[0] = c_mu * 2.0 * tke / epsilon;

            // compute derivative w.r.t. TURBULENT_KINETIC_ENERGY
            nu_t_derivatives[1] = -c_mu * std::pow(tke / epsilon, 2);
        }

        if (nu_t > mMinValue) {
            rElement.SetValue(TURBULENT_VISCOSITY, nu_t);
            rElement.SetValue(TURBULENT_VISCOSITY_DERIVATIVES, nu_t_derivatives);
        } else {
            rElement.SetValue(TURBULENT_VISCOSITY, mMinValue);
            rElement.SetValue(TURBULENT_VISCOSITY_DERIVATIVES, ZeroVector(derivatives_size));
        }
    });

    // setting wall nu_t values and their derivatives
    // TODO: for conditions
    // RansCalculationUtilities::CalculateWallTurbulentViscosity(r_model_part, mMinValue);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t derivatives for nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansNutKEpsilonDerivativesProcess::Info() const
{
    return std::string("RansNutKEpsilonDerivativesProcess");
}

void RansNutKEpsilonDerivativesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutKEpsilonDerivativesProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansNutKEpsilonDerivativesProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "min_value"       : 1e-15
        })");
    return default_parameters;
}

} // namespace Kratos.
