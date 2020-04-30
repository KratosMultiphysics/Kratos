//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah (https://github.com/sdharmin)
//                   Bence Rochlitz (https://github.com/bencerochlitz)
//
//  Supervised by:   Jordi Cotela (https://github.com/jcotela)
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/evm_k_epsilon/element_data/evm_k_epsilon_element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "evm_k_omega_element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_omega_k_element_data.h"

namespace Kratos
{
namespace EvmKOmegaElementDataUtilities
{
template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarRateVariable()
{
    return TURBULENT_KINETIC_ENERGY_RATE;
}

template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarRelaxedRateVariable()
{
    return RANS_AUXILIARY_VARIABLE_1;
}

template <unsigned int TDim>
void KElementData<TDim>::Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const NodeType& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(
            TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}
template <unsigned int TDim>
GeometryData::IntegrationMethod KElementData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateConstants(const ProcessInfo& rCurrentProcessInfo)
{
    mSigmaK = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateGaussPointData(const Vector& rShapeFunctions,
                                                 const Matrix& rShapeFunctionDerivatives,
                                                 const int Step)
{
    KRATOS_TRY

    mTurbulentKineticEnergy = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_KINETIC_ENERGY, rShapeFunctions);
    mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_VISCOSITY, rShapeFunctions);
    mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), KINEMATIC_VISCOSITY, rShapeFunctions);

    mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        this->GetGeometry(), VELOCITY, rShapeFunctionDerivatives);

    RansCalculationUtilities::CalculateGradient<TDim>(
        this->mVelocityGradient, this->GetGeometry(), VELOCITY,
        rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
array_1d<double, 3> KElementData<TDim>::CalculateEffectiveVelocity(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives) const
{
    return RansCalculationUtilities::EvaluateInPoint(this->GetGeometry(),
                                                     VELOCITY, rShapeFunctions);
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives) const
{
    return mKinematicViscosity + mSigmaK * mTurbulentKinematicViscosity;
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateReactionTerm(const Vector& rShapeFunctions,
                                                 const Matrix& rShapeFunctionDerivatives) const
{
    return std::max(mBetaStar * mTurbulentKineticEnergy / mTurbulentKinematicViscosity +
                        (2.0 / 3.0) * mVelocityDivergence,
                    0.0);
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateSourceTerm(const Vector& rShapeFunctions,
                                               const Matrix& rShapeFunctionDerivatives) const
{
    return EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);
}

// template instantiations
template class KElementData<2>;
template class KElementData<3>;

} // namespace EvmKOmegaElementDataUtilities

} // namespace Kratos