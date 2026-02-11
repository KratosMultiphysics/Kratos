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
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/serializer.h"

namespace Kratos
{

CoulombWithTensionCutOffImpl::CoulombWithTensionCutOffImpl(const Properties& rMaterialProperties)
    : mCoulombYieldSurface{rMaterialProperties}, mTensionCutOff{rMaterialProperties[GEO_TENSILE_STRENGTH]}
{
    if (rMaterialProperties.Has(GEO_ABS_YIELD_FUNCTION_TOLERANCE)) {
        mAbsoluteYieldFunctionValueTolerance = rMaterialProperties[GEO_ABS_YIELD_FUNCTION_TOLERANCE];
    }

    if (rMaterialProperties.Has(GEO_MAX_PLASTIC_ITERATIONS)) {
        mMaxNumberOfPlasticIterations = rMaterialProperties[GEO_MAX_PLASTIC_ITERATIONS];
    }
}

bool CoulombWithTensionCutOffImpl::IsAdmissibleStressState(const Geo::SigmaTau& rTrialTraction) const
{
    return IsAdmissibleStressState<>(rTrialTraction);
}

bool CoulombWithTensionCutOffImpl::IsAdmissibleStressState(const Geo::PrincipalStresses& rTrialPrincipalStresses) const
{
    return IsAdmissibleStressState<>(rTrialPrincipalStresses);
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::DoReturnMapping(const Geo::SigmaTau& rTrialTraction,
                                                            CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto sigma_tau_to_sigma_tau = [](const Geo::SigmaTau& rTraction) { return rTraction; };
    return DoReturnMapping<>(rTrialTraction, sigma_tau_to_sigma_tau, AveragingType);
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::DoReturnMapping(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                                     CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto principal_stresses_to_sigma_tau = [](const Geo::PrincipalStresses& rPrincipalStresses) {
        return StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses);
    };
    return DoReturnMapping<>(rTrialPrincipalStresses, principal_stresses_to_sigma_tau, AveragingType);
}

void CoulombWithTensionCutOffImpl::SaveKappaOfCoulombYieldSurface()
{
    mSavedKappaOfCoulombYieldSurface = mCoulombYieldSurface.GetKappa();
}

void CoulombWithTensionCutOffImpl::RestoreKappaOfCoulombYieldSurface()
{
    mCoulombYieldSurface.SetKappa(mSavedKappaOfCoulombYieldSurface);
}

template <typename StressStateType>
bool CoulombWithTensionCutOffImpl::IsAdmissibleStressState(const StressStateType& rTrialStressState) const
{
    const auto coulomb_yield_function_value = mCoulombYieldSurface.YieldFunctionValue(rTrialStressState);
    const auto tension_yield_function_value = mTensionCutOff.YieldFunctionValue(rTrialStressState);
    constexpr auto tolerance                = 1.0e-10;
    const auto     coulomb_tolerance = tolerance * (1.0 + std::abs(coulomb_yield_function_value));
    const auto     tension_tolerance = tolerance * (1.0 + std::abs(tension_yield_function_value));
    return coulomb_yield_function_value < coulomb_tolerance && tension_yield_function_value < tension_tolerance;
}

template <typename StressStateType, typename StressStateToSigmaTauFunctionType>
StressStateType CoulombWithTensionCutOffImpl::DoReturnMapping(const StressStateType& rTrialStressState,
                                                              const StressStateToSigmaTauFunctionType& rStressStateToSigmaTau,
                                                              CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto result = StressStateType{};

    const auto trial_traction = rStressStateToSigmaTau(rTrialStressState);

    auto kappa_start = mCoulombYieldSurface.GetKappa();
    for (auto counter = std::size_t{0}; counter < mMaxNumberOfPlasticIterations; ++counter) {
        if (IsStressAtTensionApexReturnZone(trial_traction)) {
            return ReturnStressAtTensionApexReturnZone(rTrialStressState);
        }

        if (IsStressAtTensionCutoffReturnZone(trial_traction)) {
            return ReturnStressAtTensionCutoffReturnZone(rTrialStressState);
        }

        if (IsStressAtCornerReturnZone(trial_traction, AveragingType)) {
            result = CalculateCornerPoint(rTrialStressState);
        } else { // Regular failure region
            result = ReturnStressAtRegularFailureZone(rTrialStressState, AveragingType);
        }

        const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrainIncrement(
                                             trial_traction, AveragingType);
        mCoulombYieldSurface.SetKappa(kappa);

        if (std::abs(mCoulombYieldSurface.YieldFunctionValue(result)) < mAbsoluteYieldFunctionValueTolerance) {
            break;
        }
    }

    return result;
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::CalculateCornerPoint(const Geo::SigmaTau&) const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    if (const auto apex = mCoulombYieldSurface.CalculateApex(); tensile_strength > apex)
        return Geo::SigmaTau{apex, 0.0};

    const auto cohesion = mCoulombYieldSurface.GetCohesion();
    const auto sin_phi  = std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto cos_phi  = std::cos(mCoulombYieldSurface.GetFrictionAngleInRadians());
    return Geo::SigmaTau{(tensile_strength - cohesion * cos_phi) / (1.0 - sin_phi),
                         (cohesion * cos_phi - tensile_strength * sin_phi) / (1.0 - sin_phi)};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::CalculateCornerPoint(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        CalculateCornerPoint(StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses)),
        rPrincipalStresses);
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionApexReturnZone(const Geo::SigmaTau& rTrialTraction) const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    return tensile_strength < mCoulombYieldSurface.CalculateApex() &&
           rTrialTraction.Sigma() - rTrialTraction.Tau() - tensile_strength > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialTraction) const
{
    const auto corner_point = CalculateCornerPoint(rTrialTraction);
    return mTensionCutOff.GetTensileStrength() < mCoulombYieldSurface.CalculateApex() &&
           corner_point.Tau() - rTrialTraction.Tau() - corner_point.Sigma() + rTrialTraction.Sigma() > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtCornerReturnZone(const Geo::SigmaTau& rTrialTraction,
                                                              CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto corner_point = CalculateCornerPoint(rTrialTraction);
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialTraction, AveragingType);
    return (rTrialTraction.Sigma() - corner_point.Sigma()) * derivative_of_flow_function[1] -
               (rTrialTraction.Tau() - corner_point.Tau()) * derivative_of_flow_function[0] >=
           0.0;
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone(const Geo::SigmaTau&) const
{
    return Geo::SigmaTau{mTensionCutOff.GetTensileStrength(), 0.0};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        ReturnStressAtTensionApexReturnZone(
            StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses)),
        rPrincipalStresses);
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTraction) const
{
    const auto derivative_of_flow_function = mTensionCutOff.DerivativeOfFlowFunction(rTraction);
    const auto lambda_tc = (mTensionCutOff.GetTensileStrength() - rTraction.Sigma() - rTraction.Tau()) /
                           (derivative_of_flow_function[0] + derivative_of_flow_function[1]);
    return rTraction + lambda_tc * derivative_of_flow_function;
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(
    const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        ReturnStressAtTensionCutoffReturnZone(
            StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses)),
        rPrincipalStresses);
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(
    const Geo::SigmaTau& rTraction, CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTraction, AveragingType);
    const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(rTraction, derivative_of_flow_function);
    return rTraction + lambda * derivative_of_flow_function;
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(
    const Geo::PrincipalStresses& rPrincipalStresses, CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        ReturnStressAtRegularFailureZone(
            StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses), AveragingType),
        rPrincipalStresses);
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
