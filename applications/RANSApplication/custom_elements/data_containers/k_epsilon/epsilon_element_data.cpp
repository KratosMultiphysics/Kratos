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
#include "custom_utilities/rans_calculation_utilities.h"
#include "element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "epsilon_element_data.h"

namespace Kratos
{
namespace KEpsilonElementData
{
template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarRateVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE_2;
}

template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarRelaxedRateVariable()
{
    return RANS_AUXILIARY_VARIABLE_2;
}

template <unsigned int TDim>
void EpsilonElementData<TDim>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void EpsilonElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mC1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    mC2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
    mCmu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mInvEpsilonSigma = 1.0 / rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
}

template <unsigned int TDim>
void EpsilonElementData<TDim>::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    double tke;

    EvaluateInPoint(this->GetGeometry(), rShapeFunctions, Step,
                    std::tie(tke, TURBULENT_KINETIC_ENERGY),
                    std::tie(mKinematicViscosity, KINEMATIC_VISCOSITY),
                    std::tie(mTurbulentKinematicViscosity, TURBULENT_VISCOSITY),
                    std::tie(mEffectiveVelocity, VELOCITY));

    mGamma = KEpsilonElementData::CalculateGamma(mCmu, tke, mTurbulentKinematicViscosity);

    mVelocityDivergence =
        GetDivergence(this->GetGeometry(), VELOCITY, rShapeFunctionDerivatives);

    CalculateGradient<TDim>(this->mVelocityGradient, this->GetGeometry(),
                            VELOCITY, rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
array_1d<double, 3> EpsilonElementData<TDim>::CalculateEffectiveVelocity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mEffectiveVelocity;
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mKinematicViscosity + mTurbulentKinematicViscosity * mInvEpsilonSigma;
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateReactionTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return std::max(mC2 * mGamma + mC1 * 2.0 * mVelocityDivergence / 3.0, 0.0);
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateSourceTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    double production = 0.0;

    production = KEpsilonElementData::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);

    production *= (mC1 * mGamma);
    return production;
}

// template instantiations

template class EpsilonElementData<2>;
template class EpsilonElementData<3>;

} // namespace KEpsilonElementData

} // namespace Kratos