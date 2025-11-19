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
            return ReturnStressAtTensionCutoffReturnZone(rTrialSigmaTau);
        }

        if (IsStressAtCornerReturnZone(rTrialSigmaTau, AveragingType)) {
            result = CalculateCornerPoint();
        } else { // Regular failure region
            result = ReturnStressAtRegularFailureZone(rTrialSigmaTau, AveragingType);
        }

        const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrain(
                                             rTrialSigmaTau, AveragingType);
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

Vector CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau) const
{
    const auto derivative_of_flow_function = mTensionCutOff.DerivativeOfFlowFunction(rSigmaTau);
    const auto lambda_tc = (mTensionCutOff.GetTensileStrength() - rSigmaTau[0] - rSigmaTau[1]) /
                           (derivative_of_flow_function[0] + derivative_of_flow_function[1]);
    return rSigmaTau + lambda_tc * derivative_of_flow_function;
}

Vector CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(const Vector& rSigmaTau,
                                                                      CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rSigmaTau, AveragingType);
    const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(rSigmaTau, derivative_of_flow_function);
    return rSigmaTau + lambda * derivative_of_flow_function;
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
