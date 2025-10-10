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

double CalculatePlasticMultiplier(const Vector& rSigmaTau,
                                  const Vector& rDerivativeOfFlowFunction,
                                  double        FrictionAngleInRadians,
                                  double Cohesion);

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
    double lambda = CalculatePlasticMultiplier(rSigmaTau, rDerivativeOfFlowFunction, FrictionAngleInRadians, Cohesion);
    return rSigmaTau + lambda * rDerivativeOfFlowFunction;
}

double CalculatePlasticMultiplier(const Vector& rSigmaTau,
                                  const Vector& rDerivativeOfFlowFunction,
                                  double        FrictionAngleInRadians,
                                  double        Cohesion)
{
    const auto sin_phi   = std::sin(FrictionAngleInRadians);
    const auto numerator = sin_phi * rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1];
    const auto lambda =
        (Cohesion * std::cos(FrictionAngleInRadians) - rSigmaTau[0] * sin_phi - rSigmaTau[1]) / numerator;
    return lambda;
}

} // namespace

namespace Kratos
{
CoulombWithTensionCutOffImpl::CoulombWithTensionCutOffImpl(double FrictionAngleInRadians,
                                                           double Cohesion,
                                                           double DilatancyAngleInRadians,
                                                           double TensileStrength)
    : mCoulombYieldSurface{FrictionAngleInRadians, Cohesion, DilatancyAngleInRadians},
      mTensionCutOff{TensileStrength}
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

Vector CoulombWithTensionCutOffImpl::DoReturnMapping(const Properties& rProperties,
                                                     const Vector&     rTrialSigmaTau,
                                                     CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    Vector result = ZeroVector(2);

    double     kappa = mEquivalentPlasticStrain;
    double tolerance = 1.0;
    while (tolerance > 1.0e-6) {
        double lambda = CalculatePlasticMultiplier(rTrialSigmaTau,
                mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType),
                mCoulombYieldSurface.GetFrictionAngleInRadians(),
                mCoulombYieldSurface.GetCohesion());
        double delta_kappa = CalculateEquivalentPlasticStrain(rTrialSigmaTau, AveragingType, lambda);
        kappa += delta_kappa;
        mEquivalentPlasticStrain = kappa;

        const auto apex =
            CalculateApex(mCoulombYieldSurface.GetFrictionAngleInRadians(), mCoulombYieldSurface.GetCohesion());

        const auto corner_point =
            CalculateCornerPoint(mCoulombYieldSurface.GetFrictionAngleInRadians(),
                                 mCoulombYieldSurface.GetCohesion(), mTensionCutOff.GetTensileStrength());

        if (IsStressAtTensionApexReturnZone(rTrialSigmaTau, mTensionCutOff.GetTensileStrength(), apex)) {
            return ReturnStressAtTensionApexReturnZone(mTensionCutOff.GetTensileStrength());
        }
        else if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau, mTensionCutOff.GetTensileStrength(), apex, corner_point)) {
            return ReturnStressAtTensionCutoffReturnZone(
                rTrialSigmaTau, mTensionCutOff.DerivativeOfFlowFunction(rTrialSigmaTau),
                mTensionCutOff.GetTensileStrength());
        }
        else if (IsStressAtCornerReturnZone(
                rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType),
                corner_point)) {
            return corner_point;
                }
        else {
            // Regular failure region
            result = ReturnStressAtRegularFailureZone(
                rTrialSigmaTau, mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType),
                mCoulombYieldSurface.GetFrictionAngleInRadians(), mCoulombYieldSurface.GetCohesion());
        }

        mCoulombYieldSurface.UpdateSurfaceProperties(rProperties[GEO_FRICTION_ANGLE],
                                                     rProperties[GEO_FRICTION_ANGLE_STRENGTH_FACTOR],
                                                     rProperties[GEO_COHESION],
                                                     rProperties[GEO_COHESION_STRENGTH_FACTOR],
                                                     rProperties[GEO_DILATANCY_ANGLE],
                                                     rProperties[GEO_DILATANCY_ANGLE_STRENGTH_FACTOR],
                                                     kappa);

        tolerance = std::abs(mCoulombYieldSurface.YieldFunctionValue(result));
    }
    return result;
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

double CoulombWithTensionCutOffImpl::CalculateEquivalentPlasticStrain(const Vector& rSigmaTau,
    CoulombYieldSurface::CoulombAveragingType AveragingType, double lambda) const
{
    Vector dGdsigma =  mCoulombYieldSurface.DerivativeOfFlowFunction(rSigmaTau, AveragingType);
    double g1 = (dGdsigma[0] + dGdsigma[1]) * 0.5;
    double g3 = (dGdsigma[0] - dGdsigma[1]) * 0.5;

    double m_mean = (g1 + g3) / 3.0;
    double m_deviatoric = std::sqrt(std::pow(g1-m_mean, 2) + std::pow(g3-m_mean, 2));
    double alpha = std::sqrt(2.0/3.0) * m_deviatoric;
    double delta_kappa = alpha * lambda;

    return delta_kappa;
}

} // namespace Kratos
