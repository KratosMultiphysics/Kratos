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
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/serializer.h"

namespace
{

using namespace Kratos;

double CalculateApex(double FrictionAngleInRadians, double Cohesion)
{
    return Cohesion / std::tan(FrictionAngleInRadians);
}

Vector CalculateCornerPoint(double FrictionAngleInRadians, double Cohesion, double TensileStrength)
{
    // Check whether the tension cut-off lies beyond the apex
    auto result = Vector{ZeroVector(2)};
    result[0]   = CalculateApex(FrictionAngleInRadians, Cohesion);
    if (TensileStrength > result[0]) return result;

    result[0] = (TensileStrength - Cohesion * std::cos(FrictionAngleInRadians)) /
                (1.0 - std::sin(FrictionAngleInRadians));
    result[1] = (Cohesion * std::cos(FrictionAngleInRadians) - TensileStrength * std::sin(FrictionAngleInRadians)) /
                (1.0 - std::sin(FrictionAngleInRadians));
    return result;
}

bool IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau, double TensileStrength, double Apex)
{
    return TensileStrength < Apex && rTrialSigmaTau[0] - rTrialSigmaTau[1] - TensileStrength > 0.0;
}

bool IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau, double TensileStrength, double Apex, const Vector& rCornerPoint)
{
    return TensileStrength < Apex &&
           rCornerPoint[1] - rTrialSigmaTau[1] - rCornerPoint[0] + rTrialSigmaTau[0] > 0.0;
}

bool IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                const Vector& rDerivativeOfFlowFunction,
                                const Vector& rCornerPoint)
{
    return (rTrialSigmaTau[0] - rCornerPoint[0]) * rDerivativeOfFlowFunction[1] -
               (rTrialSigmaTau[1] - rCornerPoint[1]) * rDerivativeOfFlowFunction[0] >=
           0.0;
}

Vector ReturnStressAtTensionApexReturnZone(double TensileStrength)
{
    auto result = Vector{ZeroVector{2}};
    result[0]   = TensileStrength;
    return result;
}

Vector ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau,
                                             const Vector& rDerivativeOfFlowFunction,
                                             double        TensileStrength)
{
    const auto lambda = (TensileStrength - rSigmaTau[0] - rSigmaTau[1]) /
                        (rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1]);
    return rSigmaTau + lambda * rDerivativeOfFlowFunction;
}

Vector ReturnStressAtRegularFailureZone(const Vector& rSigmaTau,
                                        const Vector& rDerivativeOfFlowFunction,
                                        double        FrictionAngleInRadians,
                                        double        Cohesion)
{
    const auto sin_phi   = std::sin(FrictionAngleInRadians);
    const auto numerator = sin_phi * rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1];
    const auto lambda =
        (Cohesion * std::cos(FrictionAngleInRadians) - rSigmaTau[0] * sin_phi - rSigmaTau[1]) / numerator;
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
                                                     CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto apex = CalculateApex(mCoulombYieldSurface.GetFrictionAngleInRadians(),
                                    mCoulombYieldSurface.GetCohesion());

    if (IsStressAtTensionApexReturnZone(rTrialSigmaTau, mTensionCutOff.GetTensileStrength(), apex)) {
        return ReturnStressAtTensionApexReturnZone(mTensionCutOff.GetTensileStrength());
    }

    const auto corner_point =
        CalculateCornerPoint(mCoulombYieldSurface.GetFrictionAngleInRadians(),
                             mCoulombYieldSurface.GetCohesion(), mTensionCutOff.GetTensileStrength());
    if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau, mTensionCutOff.GetTensileStrength(), apex, corner_point)) {
        return ReturnStressAtTensionCutoffReturnZone(
            rTrialSigmaTau, mTensionCutOff.DerivativeOfFlowFunction(rTrialSigmaTau),
            mTensionCutOff.GetTensileStrength());
    }

    if (IsStressAtCornerReturnZone(
            rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType),
            corner_point)) {
        return corner_point;
    }

    // Regular failure region
    return ReturnStressAtRegularFailureZone(
        rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType),
        mCoulombYieldSurface.GetFrictionAngleInRadians(), mCoulombYieldSurface.GetCohesion());
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
