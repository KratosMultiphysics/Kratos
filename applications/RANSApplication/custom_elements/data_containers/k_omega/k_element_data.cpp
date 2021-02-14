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
#include "k_element_data.h"

namespace Kratos
{
namespace KOmegaElementData
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
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_KINETIC_ENERGY_SIGMA))
        << "TURBULENT_KINETIC_ENERGY_SIGMA is not found in process info.\n";

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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_KINETIC_ENERGY))
            << "TURBULENT_KINETIC_ENERGY is not found in non-historical data "
               "value container of node with id "
            << r_node.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
            << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in non-historical data value container of node with id "
            << r_node.Id() << ".\n";

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mSigmaK = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const auto& r_properties = this->GetProperties();
    mDensity = r_properties[DENSITY];
    mKinematicViscosity = r_properties[DYNAMIC_VISCOSITY] / mDensity;
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateGaussPointData(
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

    mEffectiveKinematicViscosity = mKinematicViscosity + mSigmaK * mTurbulentKinematicViscosity;
    mReactionTerm = std::max(mBetaStar * mTurbulentKineticEnergy / mTurbulentKinematicViscosity + (2.0 / 3.0) * mVelocityDivergence, 0.0);
    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, mTurbulentKinematicViscosity);

    KRATOS_CATCH("");
}

// template instantiations
template class KElementData<2>;
template class KElementData<3>;

} // namespace KOmegaElementData

} // namespace Kratos