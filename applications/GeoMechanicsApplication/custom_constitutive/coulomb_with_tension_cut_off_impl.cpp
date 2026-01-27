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
#include "custom_constitutive/geo_principal_stresses.hpp"
#include "custom_constitutive/geo_sigma_tau.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/ublas_utilities.h"
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

bool CoulombWithTensionCutOffImpl::IsAdmissibleStressState(const Geo::SigmaTau& rTrialSigmaTau) const
{
    return IsAdmissibleStressState<>(rTrialSigmaTau);
}

bool CoulombWithTensionCutOffImpl::IsAdmissibleStressState(const Geo::PrincipalStresses& rTrialPrincipalStresses) const
{
    return IsAdmissibleStressState<>(rTrialPrincipalStresses);
}

Vector CoulombWithTensionCutOffImpl::DoReturnMapping(const Vector& rTrialSigmaTau,
                                                     CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    Vector result = ZeroVector(2);

    auto kappa_start = mCoulombYieldSurface.GetKappa();
    for (auto counter = std::size_t{0}; counter < mMaxNumberOfPlasticIterations; ++counter) {
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

        const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrainIncrement(
                                             rTrialSigmaTau, AveragingType);
        mCoulombYieldSurface.SetKappa(kappa);

        if (std::abs(mCoulombYieldSurface.YieldFunctionValue(result)) < mAbsoluteYieldFunctionValueTolerance) {
            break;
        }
    }
    return result;
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::DoReturnMapping(const Geo::SigmaTau& rTrialSigmaTau,
                                                            CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto result = Geo::SigmaTau{Vector{2, 0.0}};

    auto kappa_start = mCoulombYieldSurface.GetKappa();
    for (auto counter = std::size_t{0}; counter < mMaxNumberOfPlasticIterations; ++counter) {
        if (IsStressAtTensionApexReturnZone(rTrialSigmaTau)) {
            return ReturnStressAtTensionApexReturnZone(rTrialSigmaTau);
        }

        if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau)) {
            return ReturnStressAtTensionCutoffReturnZone(rTrialSigmaTau);
        }

        if (IsStressAtCornerReturnZone(rTrialSigmaTau, AveragingType)) {
            result = Geo::SigmaTau{CalculateCornerPoint()};
        } else { // Regular failure region
            result = ReturnStressAtRegularFailureZone(rTrialSigmaTau, AveragingType);
        }

        const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrainIncrement(
                                             rTrialSigmaTau, AveragingType);
        mCoulombYieldSurface.SetKappa(kappa);

        if (std::abs(mCoulombYieldSurface.YieldFunctionValue(result)) < mAbsoluteYieldFunctionValueTolerance) {
            break;
        }
    }

    return result;
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::DoReturnMapping(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                                     CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto result = Geo::PrincipalStresses{Vector{3, 0.0}};

    auto kappa_start = mCoulombYieldSurface.GetKappa();
    for (auto counter = std::size_t{0}; counter < mMaxNumberOfPlasticIterations; ++counter) {
        //     if (IsStressAtTensionApexReturnZone(rTrialSigmaTau)) {
        //         return ReturnStressAtTensionApexReturnZone();
        //     }
        //
        //     if (IsStressAtTensionCutoffReturnZone(rTrialSigmaTau)) {
        //         return ReturnStressAtTensionCutoffReturnZone(rTrialSigmaTau);
        //     }
        //
        //     if (IsStressAtCornerReturnZone(rTrialSigmaTau, AveragingType)) {
        //         result = CalculateCornerPoint();
        //     } else { // Regular failure region
        //         result = ReturnStressAtRegularFailureZone(rTrialSigmaTau, AveragingType);
        //     }
        //
        //     const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrainIncrement(
        //                                          rTrialSigmaTau, AveragingType);
        //     mCoulombYieldSurface.SetKappa(kappa);
        //
        //     if (std::abs(mCoulombYieldSurface.YieldFunctionValue(result)) < mAbsoluteYieldFunctionValueTolerance) {
        //         break;
        //     }
    }

    return result;
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

Vector CoulombWithTensionCutOffImpl::CalculateCornerPoint() const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    if (const auto apex = mCoulombYieldSurface.CalculateApex(); tensile_strength > apex)
        return UblasUtilities::CreateVector({apex, 0.0});

    const auto cohesion = mCoulombYieldSurface.GetCohesion();
    const auto sin_phi  = std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto cos_phi  = std::cos(mCoulombYieldSurface.GetFrictionAngleInRadians());
    return UblasUtilities::CreateVector({(tensile_strength - cohesion * cos_phi) / (1.0 - sin_phi),
                                         (cohesion * cos_phi - tensile_strength * sin_phi) / (1.0 - sin_phi)});
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
bool CoulombWithTensionCutOffImpl::IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau) const
{
    return IsStressAtTensionApexReturnZone(Geo::SigmaTau{rTrialSigmaTau});
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionApexReturnZone(const Geo::SigmaTau& rTrialSigmaTau) const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    return tensile_strength < mCoulombYieldSurface.CalculateApex() &&
           rTrialSigmaTau.sigma - rTrialSigmaTau.tau - tensile_strength > 0.0;
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
bool CoulombWithTensionCutOffImpl::IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau) const
{
    return IsStressAtTensionCutoffReturnZone(Geo::SigmaTau{rTrialSigmaTau});
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialSigmaTau) const
{
    const auto corner_point = CalculateCornerPoint();
    return mTensionCutOff.GetTensileStrength() < mCoulombYieldSurface.CalculateApex() &&
           corner_point[1] - rTrialSigmaTau.tau - corner_point[0] + rTrialSigmaTau.sigma > 0.0;
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
bool CoulombWithTensionCutOffImpl::IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                                              CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    return IsStressAtCornerReturnZone(Geo::SigmaTau{rTrialSigmaTau}, AveragingType);
}

bool CoulombWithTensionCutOffImpl::IsStressAtCornerReturnZone(const Geo::SigmaTau& rTrialSigmaTau,
                                                              CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto corner_point = CalculateCornerPoint();
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialSigmaTau, AveragingType);
    return (rTrialSigmaTau.sigma - corner_point[0]) * derivative_of_flow_function[1] -
               (rTrialSigmaTau.tau - corner_point[1]) * derivative_of_flow_function[0] >=
           0.0;
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
Vector CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone() const
{
    const auto unused_sigma_tau = Geo::SigmaTau{Vector{ZeroVector{2}}};
    return ReturnStressAtTensionApexReturnZone(unused_sigma_tau).CopyTo<Vector>();
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone(const Geo::SigmaTau&) const
{
    return Geo::SigmaTau{UblasUtilities::CreateVector({mTensionCutOff.GetTensileStrength(), 0.0})};
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
Vector CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau) const
{
    return ReturnStressAtTensionCutoffReturnZone(Geo::SigmaTau{rSigmaTau}).CopyTo<Vector>();
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rSigmaTau) const
{
    const auto derivative_of_flow_function = mTensionCutOff.DerivativeOfFlowFunction(rSigmaTau);
    const auto lambda_tc = (mTensionCutOff.GetTensileStrength() - rSigmaTau.sigma - rSigmaTau.tau) /
                           (derivative_of_flow_function[0] + derivative_of_flow_function[1]);
    return Geo::SigmaTau{Vector{rSigmaTau.CopyTo<Vector>() + lambda_tc * derivative_of_flow_function}};
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
Vector CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(const Vector& rSigmaTau,
                                                                      CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    return ReturnStressAtRegularFailureZone(Geo::SigmaTau{rSigmaTau}, AveragingType).CopyTo<Vector>();
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(
    const Geo::SigmaTau& rSigmaTau, CoulombYieldSurface::CoulombAveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rSigmaTau, AveragingType);
    const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(rSigmaTau, derivative_of_flow_function);
    return Geo::SigmaTau{Vector{rSigmaTau.CopyTo<Vector>() + lambda * derivative_of_flow_function}};
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
