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

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_element_data.h"

namespace Kratos
{
namespace KOmegaElementData
{
template <unsigned int TDim>
const Variable<double>& OmegaElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_BETA))
        << "TURBULENCE_RANS_BETA is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_GAMMA))
        << "TURBULENCE_RANS_GAMMA is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);

        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_KINETIC_ENERGY))
            << "TURBULENT_KINETIC_ENERGY is not found in non-historical data "
               "value container of node with id "
            << r_node.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
            << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in non-historical data value container of node with id "
            << r_node.Id() << ".\n";

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mBeta = rCurrentProcessInfo[TURBULENCE_RANS_BETA];
    mGamma = rCurrentProcessInfo[TURBULENCE_RANS_GAMMA];
    mSigmaOmega = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA];

    const auto& r_properties = this->GetProperties();
    mDensity = r_properties[DENSITY];
    mKinematicViscosity = r_properties[DYNAMIC_VISCOSITY] / mDensity;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mEffectiveVelocity, VELOCITY));

    FluidCalculationUtilities::EvaluateGradientInPoint(
        this->GetGeometry(), rShapeFunctionDerivatives,
        std::tie(mVelocityGradient, VELOCITY));

    double tke_old, omega_old;
    FluidCalculationUtilities::EvaluateNonHistoricalInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(tke_old, TURBULENT_KINETIC_ENERGY),
        std::tie(omega_old, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    mTurbulentKinematicViscosity = 1e-12;
    if (tke_old > 0.0 && omega_old > 0.0) {
        mTurbulentKinematicViscosity = std::max(tke_old / omega_old, 1e-12);
    }

    mVelocityDivergence = CalculateMatrixTrace<TDim>(mVelocityGradient);

    mEffectiveKinematicViscosity = mKinematicViscosity + mTurbulentKinematicViscosity * mSigmaOmega;
    mReactionTerm = std::max(mBeta * mTurbulentKineticEnergy / mTurbulentKinematicViscosity + mGamma * 2.0 * mVelocityDivergence / 3.0, 0.0);
    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, 1.0) *  (mGamma);

    KRATOS_CATCH("");
}

// template instantiations

template class OmegaElementData<2>;
template class OmegaElementData<3>;

} // namespace KOmegaElementData

} // namespace Kratos