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
#include <functional>
#include <limits>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element_data/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_k_omega_sst_update_process.h"

namespace Kratos
{
RansNutKOmegaSSTUpdateProcess::RansNutKOmegaSSTUpdateProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "a1"              : 0.31,
            "b1"              : 1.0,
            "beta_star"       : 0.09,
            "min_value"       : 1e-15
        })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mMinValue = rParameters["min_value"].GetDouble();
    mA1 = rParameters["a1"].GetDouble();
    mB1 = rParameters["b1"].GetDouble();
    mBetaStar = rParameters["beta_star"].GetDouble();

    KRATOS_CATCH("");
}

RansNutKOmegaSSTUpdateProcess::RansNutKOmegaSSTUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double A1,
    const double BetaStar,
    const double MinValue,
    const int EchoLevel)
: mrModel(rModel),
  mModelPartName(rModelPartName),
  mA1(A1),
  mBetaStar(BetaStar),
  mMinValue(MinValue),
  mEchoLevel(EchoLevel)
{
}

int RansNutKOmegaSSTUpdateProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);

    return 0;

    KRATOS_CATCH("");
}

void RansNutKOmegaSSTUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized) {
        this->Execute();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansNutKOmegaSSTUpdateProcess::ExecuteInitialize()
{
    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(NUMBER_OF_NEIGHBOUR_ELEMENTS,
                                                   r_model_part.Nodes());

    const int number_of_elements = r_model_part.NumberOfElements();
#pragma omp parallel for
    for (int i_elem = 0; i_elem < number_of_elements; ++i_elem) {
        auto& r_element = *(r_model_part.ElementsBegin() + i_elem);
        auto& r_geometry = r_element.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            auto& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) += 1;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour elements in " << mModelPartName << ".\n";
}

void RansNutKOmegaSSTUpdateProcess::Execute()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    auto& r_nodes = r_model_part.Nodes();
    VariableUtils().SetHistoricalVariableToZero(TURBULENT_VISCOSITY, r_nodes);

    auto& r_elements = r_model_part.Elements();
    const int number_of_elements = r_elements.size();

    std::function<double(const Element&)> nut_calculation_method;
    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    if (domain_size == 2) {
        nut_calculation_method = [this](const Element& rElement) {
            return this->CalculateElementNuT<2>(rElement);
        };
    } else if (domain_size == 3) {
        nut_calculation_method = [this](const Element& rElement) {
            return this->CalculateElementNuT<3>(rElement);
        };
    } else {
        KRATOS_ERROR << "Unsupported domain size.";
    }

#pragma omp parallel for
    for (int i_elem = 0; i_elem < number_of_elements; ++i_elem) {
        auto& r_element = *(r_elements.begin() + i_elem);
        const double nut = nut_calculation_method(r_element);

        auto& r_geometry = r_element.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            auto& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) += nut;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);

    const int number_of_nodes = r_nodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        auto& r_node = *(r_nodes.begin() + i_node);
        const double number_of_neighbour_elements =
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        double& nut = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        nut = std::max(nut / number_of_neighbour_elements, mMinValue);
        double& nu = r_node.FastGetSolutionStepValue(VISCOSITY);
        nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) + nut;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t for nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double RansNutKOmegaSSTUpdateProcess::CalculateElementNuT(
    const Element& rElement) const
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    const GeometryType& r_geometry = rElement.GetGeometry();

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    GeometryType::ShapeFunctionsGradientsType shape_derivatives;
    CalculateGeometryData(r_geometry, rElement.GetIntegrationMethod(),
                          gauss_weights, shape_functions, shape_derivatives);
    const int num_gauss_points = gauss_weights.size();

    BoundedMatrix<double, TDim, TDim> velocity_gradient;

    double nut = 0.0;

    for (int g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& r_gauss_shape_functions = row(shape_functions, g);

        const double tke = EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY,
                                           r_gauss_shape_functions);
        const double omega = EvaluateInPoint(
            r_geometry, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_gauss_shape_functions);
        const double nu =
            EvaluateInPoint(r_geometry, KINEMATIC_VISCOSITY, r_gauss_shape_functions);
        const double y = EvaluateInPoint(r_geometry, DISTANCE, r_gauss_shape_functions);

        CalculateGradient<TDim>(velocity_gradient, r_geometry, VELOCITY, r_shape_derivatives);

        const double f_2 = KOmegaSSTElementData::CalculateF2(tke, omega, nu, y, mBetaStar);

        const BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient =
            (velocity_gradient + trans(velocity_gradient)) * 0.5;

        const double t = norm_frobenius(symmetric_velocity_gradient) * 1.414;

        nut += KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(
            tke, omega, t, f_2, mA1);
    }

    nut /= static_cast<double>(num_gauss_points);

    return nut;

    KRATOS_CATCH("");
}

std::string RansNutKOmegaSSTUpdateProcess::Info() const
{
    return std::string("RansNutKOmegaSSTUpdateProcess");
}

void RansNutKOmegaSSTUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutKOmegaSSTUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

// template instantiations
template double RansNutKOmegaSSTUpdateProcess::CalculateElementNuT<2>(const Element&) const;
template double RansNutKOmegaSSTUpdateProcess::CalculateElementNuT<3>(const Element&) const;

} // namespace Kratos.
