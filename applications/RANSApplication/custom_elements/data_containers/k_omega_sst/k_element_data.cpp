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

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_element_data.h"

namespace Kratos
{
namespace KOmegaSSTElementData
{
template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim>
void KElementData<TDim>::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA_1))
        << "TURBULENT_KINETIC_ENERGY_SIGMA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA_2))
        << "TURBULENT_KINETIC_ENERGY_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (const auto& r_node : r_geometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void KElementData<TDim>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput = CalculateEffectiveViscosity(rCurrentProcessInfo);
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateEffectiveViscosity(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    // only needs TURBULENT_VISCOSITY, so no checks were done since this will be
    // only an internal call

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    typename GeometryType::ShapeFunctionsGradientsType shape_derivatives;
    CalculateGeometryData(this->GetGeometry(), GeometryData::IntegrationMethod::GI_GAUSS_1,
                          gauss_weights, shape_functions, shape_derivatives);
    const int num_gauss_points = gauss_weights.size();

    BoundedMatrix<double, TDim, TDim> velocity_gradient;

    double tke, omega, nu, y;
    const double rho = this->GetProperties().GetValue(DENSITY);
    const double a1 =  rCurrentProcessInfo[TURBULENCE_RANS_A1];
    const double beta_star = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    double nu_t = 0.0;

    for (int g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& r_gauss_shape_functions = row(shape_functions, g);

        auto& cl_parameters = this->GetConstitutiveLawParameters();
        cl_parameters.SetShapeFunctionsValues(r_gauss_shape_functions);

        this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
        nu /= rho;

        FluidCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), r_gauss_shape_functions,
            std::tie(tke, TURBULENT_KINETIC_ENERGY),
            std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
            std::tie(y, DISTANCE)
        );

        FluidCalculationUtilities::EvaluateGradientInPoint(
            this->GetGeometry(), r_shape_derivatives,
            std::tie(velocity_gradient, VELOCITY));

        const double f_2 = KOmegaSSTElementData::CalculateF2(tke, omega, nu, y, beta_star);

        const BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient =
            (velocity_gradient + trans(velocity_gradient)) * 0.5;

        const double t = norm_frobenius(symmetric_velocity_gradient) * 1.414;

        nu_t += KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(
            tke, omega, t, f_2, a1);
    }

    return nu_t / num_gauss_points;

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mSigmaK1 = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA_1];
    mSigmaK2 = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA_2];
    mSigmaOmega2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mDensity = this->GetProperties().GetValue(DENSITY);
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(rShapeFunctions);

    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, mKinematicViscosity);
    mKinematicViscosity /= mDensity;

    const auto& r_geometry = this->GetGeometry();

    FluidCalculationUtilities::EvaluateInPoint(
        r_geometry, rShapeFunctions, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRate, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mTurbulentKinematicViscosity, TURBULENT_VISCOSITY),
        std::tie(mWallDistance, DISTANCE),
        std::tie(mEffectiveVelocity, VELOCITY));

    KRATOS_ERROR_IF(mWallDistance < 0.0) << "Wall distance is negative at " << r_geometry;

    FluidCalculationUtilities::EvaluateGradientInPoint(
        this->GetGeometry(), rShapeFunctionDerivatives,
        std::tie(mTurbulentKineticEnergyGradient, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRateGradient, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mVelocityGradient, VELOCITY));

    mCrossDiffusion = KOmegaSSTElementData::CalculateCrossDiffusionTerm<TDim>(
        mSigmaOmega2, mTurbulentSpecificEnergyDissipationRate,
        mTurbulentKineticEnergyGradient, mTurbulentSpecificEnergyDissipationRateGradient);

    const double f_1 = KOmegaSSTElementData::CalculateF1(
        mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate,
        mKinematicViscosity, mWallDistance, mBetaStar, mCrossDiffusion, mSigmaOmega2);

    mBlendedSimgaK =
        KOmegaSSTElementData::CalculateBlendedPhi(mSigmaK1, mSigmaK2, f_1);

    mVelocityDivergence = CalculateMatrixTrace<TDim>(mVelocityGradient);

    mEffectiveKinematicViscosity = mKinematicViscosity + mBlendedSimgaK * mTurbulentKinematicViscosity;
    mReactionTerm = std::max(mBetaStar * mTurbulentKineticEnergy / mTurbulentKinematicViscosity + (2.0 / 3.0) * mVelocityDivergence, 0.0);
    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, mTurbulentKinematicViscosity);

    KRATOS_CATCH("");
}

// template instantiations
template class KElementData<2>;
template class KElementData<3>;

} // namespace KOmegaSSTElementData

} // namespace Kratos