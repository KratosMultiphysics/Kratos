// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_constitutive/coulomb_with_tension_cut_off_impl.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/serializer.h"

namespace
{

using namespace Kratos;

Vector ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau,
                                             const Vector& rDerivativeOfFlowFunction,
                                             double        TensileStrength)
{
    const auto lambda = (TensileStrength - rSigmaTau[0] - rSigmaTau[1]) /
                        (rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1]);
    return rSigmaTau + lambda * rDerivativeOfFlowFunction;
}

Vector ReturnStressAtRegularFailureZone(const Vector&        rSigmaTau,
                                        CoulombYieldSurface& rCoulombYieldSurface,
                                        const Vector&        rDerivativeOfFlowFunction)
{
    const auto lambda = rCoulombYieldSurface.CalculatePlasticMultiplier(rSigmaTau, rDerivativeOfFlowFunction);
    return rSigmaTau + lambda * rDerivativeOfFlowFunction;
}

} // namespace

namespace Kratos
{

CoulombWithTensionCutOffImpl::CoulombWithTensionCutOffImpl(const Properties& rMaterialProperties)
    : mCoulombYieldSurface{rMaterialProperties}, mTensionCutOff{rMaterialProperties[GEO_TENSILE_STRENGTH]}
{
}

bool CoulombWithTensionCutOffImpl::IsAdmissibleSigmaTau(const Vector& rTrialSigmaTau) const
{
    const auto coulomb_yield_function_value = mCoulombYieldSurface.YieldFunctionValue(rTrialSigmaTau);
    const auto     tension_yield_function_value = mTensionCutOff.YieldFunctionValue(rTrialSigmaTau);
    constexpr auto tolerance                    = 1.0e-10;
    const auto     coulomb_tolerance = tolerance * (1.0 + std::abs(coulomb_yield_function_value));
    const auto     tension_tolerance = tolerance * (1.0 + std::abs(tension_yield_function_value));
    return coulomb_yield_function_value < coulomb_tolerance && tension_yield_function_value < tension_tolerance;
}

Vector CoulombWithTensionCutOffImpl::DoReturnMapping(const Vector& rTrialSigmaTau,
                                                     CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    Vector     result      = ZeroVector(2);
    const auto kappa_start = mCoulombYieldSurface.GetKappa();

    for (unsigned int counter = 0; counter < mCoulombYieldSurface.GetMaxIterations(); ++counter) {
        if (IsStressAtTensionApexReturnZone(rTrialSigmaTau)) {
            return ReturnStressAtTensionApexReturnZone();
        }

        if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau)) {
            return ReturnStressAtTensionCutoffReturnZone(
                rTrialSigmaTau, mTensionCutOff.DerivativeOfFlowFunction(rTrialSigmaTau),
                mTensionCutOff.GetTensileStrength());
        }

        if (IsStressAtCornerReturnZone(rTrialSigmaTau, AveragingType)) {
            result = CalculateCornerPoint();
        } else { // Regular failure region
            result = ReturnStressAtRegularFailureZone(
                rTrialSigmaTau, mCoulombYieldSurface,
                mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType));
        }

        const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(
            rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType));
        const double kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrain(
                                               rTrialSigmaTau, AveragingType, lambda);
        mCoulombYieldSurface.SetKappa(kappa);

        double error = std::abs(mCoulombYieldSurface.YieldFunctionValue(result));
        if (error < mCoulombYieldSurface.GetConvergenceTolerance()) break;
    }
    return result;
}

Vector CoulombWithTensionCutOffImpl::CalculateCornerPoint() const
{
    // Check whether the tension cut-off lies beyond the apex
    auto result                 = Vector{ZeroVector(2)};
    result[0]                   = mCoulombYieldSurface.CalculateApex();
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    if (tensile_strength > result[0]) return result;

    const auto cohesion = mCoulombYieldSurface.GetCohesion();
    const auto sin_phi  = std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto cos_phi  = std::cos(mCoulombYieldSurface.GetFrictionAngleInRadians());
    result[0]           = (tensile_strength - cohesion * cos_phi) / (1.0 - sin_phi);
    result[1]           = (cohesion * cos_phi - tensile_strength * sin_phi) / (1.0 - sin_phi);
    return result;
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau) const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    return tensile_strength < mCoulombYieldSurface.CalculateApex() &&
           rTrialSigmaTau[0] - rTrialSigmaTau[1] - tensile_strength > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau) const
{
    const auto corner_point = CalculateCornerPoint();
    return mTensionCutOff.GetTensileStrength() < mCoulombYieldSurface.CalculateApex() &&
           corner_point[1] - rTrialSigmaTau[1] - corner_point[0] + rTrialSigmaTau[0] > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                                              CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto corner_point = CalculateCornerPoint();
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType);
    return (rTrialSigmaTau[0] - corner_point[0]) * derivative_of_flow_function[1] -
               (rTrialSigmaTau[1] - corner_point[1]) * derivative_of_flow_function[0] >=
           0.0;
}

Vector CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone() const
{
    return UblasUtilities::CreateVector({mTensionCutOff.GetTensileStrength(), 0.0});
}

void CoulombWithTensionCutOffImpl::save(Serializer& rSerializer) const
{
    rSerializer.save("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.save("TensionCutOff", mTensionCutOff);
}

void CoulombWithTensionCutOffImpl::load(Serializer& rSerializer)
{
    rSerializer.load("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.load("TensionCutOff", mTensionCutOff);
}

} // namespace Kratos
