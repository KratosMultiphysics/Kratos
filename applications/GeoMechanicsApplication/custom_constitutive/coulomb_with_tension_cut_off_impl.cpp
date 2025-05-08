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

double CalculateApex(double FrictionAngle, double Cohesion)
{
    return Cohesion / std::tan(FrictionAngle);
}

Vector CalculateCornerPoint(double FrictionAngle, double Cohesion, double TensileStrength)
{
    // Check whether the tension cut-off lies beyond the apex
    auto result = Vector{ZeroVector(2)};
    result[0]   = CalculateApex(FrictionAngle, Cohesion);
    if (TensileStrength > result[0]) return result;

    result[0] = (TensileStrength - Cohesion * std::cos(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    result[1] = (Cohesion * std::cos(FrictionAngle) - TensileStrength * std::sin(FrictionAngle)) /
                (1.0 - std::sin(FrictionAngle));
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

bool IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau, double DilatancyAngle, const Vector& rCornerPoint)
{
    return rTrialSigmaTau[0] - rCornerPoint[0] - (rTrialSigmaTau[1] - rCornerPoint[1]) * std::sin(DilatancyAngle) >= 0.0;
}

Vector ReturnStressAtTensionApexReturnZone(double TensileStrength)
{
    Vector result{2};
    result[0] = TensileStrength;
    result[1] = 0.0;
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
    const auto cof1      = std::sin(FrictionAngleInRadians);
    const auto cof2      = Cohesion * std::cos(FrictionAngleInRadians);
    const auto numerator = cof1 * rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1];
    const auto lambda    = (cof2 - rSigmaTau[0] * cof1 - rSigmaTau[1]) / numerator;
    return rSigmaTau + lambda * rDerivativeOfFlowFunction;
}

} // namespace

namespace Kratos
{
CoulombWithTensionCutOffImpl::CoulombWithTensionCutOffImpl(double FrictionAngleInRad,
                                                           double Cohesion,
                                                           double DilatationAngleInRad,
                                                           double TensileStrength)
    : mCoulombYieldSurface{FrictionAngleInRad, Cohesion, DilatationAngleInRad}, mTensionCutOff{TensileStrength}
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

Vector CoulombWithTensionCutOffImpl::DoReturnMapping(const Properties& rProperties, const Vector& rTrialSigmaTau) const
{
    Vector result;
    const auto apex = CalculateApex(ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties),
                                    ConstitutiveLawUtilities::GetCohesion(rProperties));
    const auto corner_point = CalculateCornerPoint(
        ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties),
        ConstitutiveLawUtilities::GetCohesion(rProperties), rProperties[GEO_TENSILE_STRENGTH]);

    if (IsStressAtTensionApexReturnZone(rTrialSigmaTau, rProperties[GEO_TENSILE_STRENGTH], apex)) {
        return ReturnStressAtTensionApexReturnZone(rProperties[GEO_TENSILE_STRENGTH]);
    }

    if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau, rProperties[GEO_TENSILE_STRENGTH], apex, corner_point)) {
        return ReturnStressAtTensionCutoffReturnZone(
            rTrialSigmaTau, mTensionCutOff.DerivativeOfFlowFunction(rTrialSigmaTau),
            rProperties[GEO_TENSILE_STRENGTH]);
    }

    if (IsStressAtCornerReturnZone(
            rTrialSigmaTau, MathUtils<>::DegreesToRadians(rProperties[GEO_DILATANCY_ANGLE]), corner_point)) {
        return corner_point;
    }

    // Regular failure region
    return ReturnStressAtRegularFailureZone(
        rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau),
        ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties),
        ConstitutiveLawUtilities::GetCohesion(rProperties));
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
