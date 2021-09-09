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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_element_data.h"

namespace Kratos
{
namespace KEpsilonElementData
{
template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim>
void KElementData<TDim>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mCmu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mInvTkeSigma = 1.0 / rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    mDensity = this->GetProperties().GetValue(DENSITY);

    KRATOS_CATCH("");
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

    double tke;
    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions, Step,
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
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
array_1d<double, 3> KElementData<TDim>::CalculateEffectiveVelocity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mEffectiveVelocity;
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mKinematicViscosity + mTurbulentKinematicViscosity * mInvTkeSigma;
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateReactionTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return std::max(mGamma + (2.0 / 3.0) * mVelocityDivergence, 0.0);
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateSourceTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    double production = 0.0;

    production = KEpsilonElementData::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);

    return production;
}

// template instantiations
template class KElementData<2>;
template class KElementData<3>;

} // namespace KEpsilonElementData

} // namespace Kratos