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
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"

namespace
{

using namespace Kratos;

template <typename YieldSurfaceType>
auto CalculatePrincipalStressCorrection(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                        const Geo::PrincipalStresses::AveragingType AveragingType,
                                        const Matrix&           rElasticConstitutiveTensor,
                                        const YieldSurfaceType& rYieldSurface)
{
    const auto dG_dSigma = rYieldSurface.DerivativeOfFlowFunction(rTrialPrincipalStresses, AveragingType);
    return Geo::PrincipalStresses{prod(subrange(rElasticConstitutiveTensor, 0, 3, 0, 3), dG_dSigma)};
}

} // namespace

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
                                                            const Matrix& rElasticConstitutiveTensor,
                                                            Geo::PrincipalStresses::AveragingType AveragingType)
{
    auto sigma_tau_to_sigma_tau = [](const Geo::SigmaTau& rTraction) { return rTraction; };
    return DoReturnMapping<>(rTrialTraction, sigma_tau_to_sigma_tau, rElasticConstitutiveTensor, AveragingType);
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::DoReturnMapping(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                                     const Matrix& rElasticConstitutiveTensor,
                                                                     Geo::PrincipalStresses::AveragingType AveragingType)
{
    auto principal_stresses_to_sigma_tau = [](const Geo::PrincipalStresses& rPrincipalStresses) {
        return StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses);
    };
    return DoReturnMapping<>(rTrialPrincipalStresses, principal_stresses_to_sigma_tau,
                             rElasticConstitutiveTensor, AveragingType);
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
                                                              const Matrix& rElasticConstitutiveTensor,
                                                              Geo::PrincipalStresses::AveragingType AveragingType)
{
    auto result = StressStateType{};

    const auto trial_traction = rStressStateToSigmaTau(rTrialStressState);

    auto kappa_start = mCoulombYieldSurface.GetKappa();
    for (auto counter = std::size_t{0}; counter < mMaxNumberOfPlasticIterations; ++counter) {
        if (IsStressAtTensionApexReturnZone(trial_traction)) {
            return ReturnStressAtTensionApexReturnZone(rTrialStressState);
        }

        if (IsStressAtTensionCutoffReturnZone(trial_traction)) {
            return ReturnStressAtTensionCutoffReturnZone(rTrialStressState,
                                                         rElasticConstitutiveTensor, AveragingType);
        }

        //if(

        if (IsStressAtCornerReturnZone(trial_traction, AveragingType)) {
            result = ReturnStressAtCornerPoint(rTrialStressState, rElasticConstitutiveTensor, AveragingType);
        } else { // Regular failure region
            result = ReturnStressAtRegularFailureZone(rTrialStressState, rElasticConstitutiveTensor, AveragingType);
        }

        const auto kappa = kappa_start + mCoulombYieldSurface.CalculateEquivalentPlasticStrainIncrement(
                                             rTrialStressState, rElasticConstitutiveTensor, AveragingType);
        mCoulombYieldSurface.SetKappa(kappa);

        if (std::abs(mCoulombYieldSurface.YieldFunctionValue(result)) < mAbsoluteYieldFunctionValueTolerance) {
            break;
        }
    }

    return result;
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::CalculateCornerPoint() const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    if (const auto apex = mCoulombYieldSurface.CalculateApex(); tensile_strength > apex.Sigma())
        return apex;

    const auto cohesion = mCoulombYieldSurface.GetCohesion();
    const auto sin_phi  = std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto cos_phi  = std::cos(mCoulombYieldSurface.GetFrictionAngleInRadians());
    return Geo::SigmaTau{(tensile_strength - cohesion * cos_phi) / (1.0 - sin_phi),
                         (cohesion * cos_phi - tensile_strength * sin_phi) / (1.0 - sin_phi)};
}

Geo::PQ CoulombWithTensionCutOffImpl::CalculateCapCornerPoint() const
{
    auto       result = Vector(2);
    const auto b1     = 1.0 / std::pow(mCompressionCapYieldSurface->GetCapSize(), 2);
    const auto c1     = std::pow(mCompressionCapYieldSurface->GetPreconsolidationStress(), 2);
    const auto sin_phi = std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto cos_phi = std::cos(mCoulombYieldSurface.GetFrictionAngleInRadians());
    const auto a2 = 6.0 * sin_phi / (3.0 - sin_phi);
    const auto c2 =
        -6.0 * mCoulombYieldSurface.GetCohesion() * cos_phi / (3.0 - sin_phi);
    const auto A     = 1.0 + b1 * a2 * a2;
    const auto B     = -2.0 * b1 * a2 * c2;
    const auto C     = b1 * c2 * c2 - c1;
    const auto delta = B * B - 4.0 * A * C;
    if (delta > 0.0) {
        const auto p1 = (-B + std::sqrt(delta)) / (2.0 * A);
        const auto p2 = (-B - std::sqrt(delta)) / (2.0 * A);
        result[0]     = std::min(p1, p2);
        result[1]     = a2 * result[0] - c2;
        return Geo::PQ{result};
    }
    KRATOS_ERROR << "Failed to calculate the cap corner point. No intersection with Coulomb yield "
                    "surface was found.\n ";
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionApexReturnZone(const Geo::SigmaTau& rTrialTraction) const
{
    const auto tensile_strength = mTensionCutOff.GetTensileStrength();
    return tensile_strength < mCoulombYieldSurface.CalculateApex().Sigma() &&
           rTrialTraction.Sigma() - rTrialTraction.Tau() - tensile_strength > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialTraction) const
{
    const auto corner_point = CalculateCornerPoint();
    return mTensionCutOff.GetTensileStrength() < mCoulombYieldSurface.CalculateApex().Sigma() &&
           corner_point.Tau() - rTrialTraction.Tau() - corner_point.Sigma() + rTrialTraction.Sigma() > 0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtCornerReturnZone(const Geo::SigmaTau& rTrialTraction,
                                                              Geo::PrincipalStresses::AveragingType AveragingType) const
{
    const auto corner_point = CalculateCornerPoint();
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialTraction, AveragingType);
    return (rTrialTraction.Sigma() - corner_point.Sigma()) * derivative_of_flow_function[1] -
               (rTrialTraction.Tau() - corner_point.Tau()) * derivative_of_flow_function[0] >=
           0.0;
}

bool CoulombWithTensionCutOffImpl::IsStressAtCompressionCapReturnZone(const Geo::PQ& rTrialPQ) const
{
    const auto cap_corner_point = CalculateCapCornerPoint();
    const auto cof = cap_corner_point.P() * std::pow(mCompressionCapYieldSurface->GetCapSize(), 2) / cap_corner_point.Q();
    return mCompressionCapYieldSurface->YieldFunctionValue(rTrialPQ) > 0.0 &&
    ((rTrialPQ.Q() - cap_corner_point.Q()) - cof * (rTrialPQ.P() - cap_corner_point.P())) < 0.0;
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone(const Geo::SigmaTau&) const
{
    return Geo::SigmaTau{mTensionCutOff.GetTensileStrength(), 0.0};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtTensionApexReturnZone(
    const Geo::PrincipalStresses& rTrialPrincipalStresses) const
{
    return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        ReturnStressAtTensionApexReturnZone(
            StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rTrialPrincipalStresses)),
        rTrialPrincipalStresses);
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(
    const Geo::SigmaTau&                  rTrialTraction,
    const Matrix&                         rElasticConstitutiveTensor,
    Geo::PrincipalStresses::AveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mTensionCutOff.DerivativeOfFlowFunction(rTrialTraction, AveragingType);
    const auto lambda = mTensionCutOff.CalculatePlasticMultiplier(
        rTrialTraction, derivative_of_flow_function, rElasticConstitutiveTensor);
    return rTrialTraction +
           Geo::SigmaTau{lambda * prod(rElasticConstitutiveTensor, derivative_of_flow_function)};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtTensionCutoffReturnZone(
    const Geo::PrincipalStresses&         rTrialPrincipalStresses,
    const Matrix&                         rElasticConstitutiveTensor,
    Geo::PrincipalStresses::AveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mTensionCutOff.DerivativeOfFlowFunction(rTrialPrincipalStresses, AveragingType);
    const auto lambda = mTensionCutOff.CalculatePlasticMultiplier(
        rTrialPrincipalStresses, derivative_of_flow_function, rElasticConstitutiveTensor);
    return rTrialPrincipalStresses +
           Geo::PrincipalStresses{lambda * prod(subrange(rElasticConstitutiveTensor, 0, 3, 0, 3),
                                                derivative_of_flow_function)};
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(
    const Geo::SigmaTau&                  rTrialTraction,
    const Matrix&                         rElasticConstitutiveTensor,
    Geo::PrincipalStresses::AveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialTraction, AveragingType);
    const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(
        rTrialTraction, derivative_of_flow_function, rElasticConstitutiveTensor);
    return rTrialTraction +
           Geo::SigmaTau{lambda * prod(rElasticConstitutiveTensor, derivative_of_flow_function)};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtRegularFailureZone(
    const Geo::PrincipalStresses&         rTrialPrincipalStresses,
    const Matrix&                         rElasticConstitutiveTensor,
    Geo::PrincipalStresses::AveragingType AveragingType) const
{
    const auto derivative_of_flow_function =
        mCoulombYieldSurface.DerivativeOfFlowFunction(rTrialPrincipalStresses, AveragingType);
    const auto lambda = mCoulombYieldSurface.CalculatePlasticMultiplier(
        rTrialPrincipalStresses, derivative_of_flow_function, rElasticConstitutiveTensor);
    return rTrialPrincipalStresses +
           Geo::PrincipalStresses{lambda * prod(subrange(rElasticConstitutiveTensor, 0, 3, 0, 3),
                                                derivative_of_flow_function)};
}

Geo::PrincipalStresses CoulombWithTensionCutOffImpl::ReturnStressAtCornerPoint(
    const Geo::PrincipalStresses&         rTrialPrincipalStresses,
    const Matrix&                         rElasticConstitutiveTensor,
    Geo::PrincipalStresses::AveragingType AveragingType) const
{
    if (const auto apex = mCoulombYieldSurface.CalculateApex();
        mTensionCutOff.GetTensileStrength() > apex.Sigma())
        return StressStrainUtilities::TransformSigmaTauToPrincipalStresses(apex, rTrialPrincipalStresses);

    const auto principal_stress_correction_Coulomb = CalculatePrincipalStressCorrection(
        rTrialPrincipalStresses, AveragingType, rElasticConstitutiveTensor, mCoulombYieldSurface);
    const auto traction_correction_Coulomb =
        StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress_correction_Coulomb);

    const auto principal_stress_correction_tension_cut_off = CalculatePrincipalStressCorrection(
        rTrialPrincipalStresses, AveragingType, rElasticConstitutiveTensor, mTensionCutOff);
    const auto traction_correction_tension_cut_off =
        StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_stress_correction_tension_cut_off);

    const auto v = std::array{std::sin(mCoulombYieldSurface.GetFrictionAngleInRadians()), 1.0};
    const auto A =
        UblasUtilities::CreateMatrix({{MathUtils<>::Dot(traction_correction_Coulomb.Values(), v),
                                       MathUtils<>::Dot(traction_correction_tension_cut_off.Values(), v)},
                                      {principal_stress_correction_Coulomb.Values()[0],
                                       principal_stress_correction_tension_cut_off.Values()[0]}});

    auto A_inverse   = Matrix(2, 2);
    auto determinant = 0.0;
    MathUtils<>::InvertMatrix2(A, A_inverse, determinant);

    const auto b = UblasUtilities::CreateVector(
        {-1.0 * mCoulombYieldSurface.YieldFunctionValue(rTrialPrincipalStresses),
         -1.0 * mTensionCutOff.YieldFunctionValue(rTrialPrincipalStresses)});
    const auto plastic_multipliers = Vector{prod(A_inverse, b)};

    return rTrialPrincipalStresses +
           Geo::PrincipalStresses{
               plastic_multipliers[0] * principal_stress_correction_Coulomb.Values() +
               plastic_multipliers[1] * principal_stress_correction_tension_cut_off.Values()};
}

Geo::SigmaTau CoulombWithTensionCutOffImpl::ReturnStressAtCornerPoint(
    const Geo::SigmaTau&, const Matrix&, Geo::PrincipalStresses::AveragingType AveragingType) const
{
    KRATOS_DEBUG_ERROR_IF(AveragingType != Geo::PrincipalStresses::AveragingType::NO_AVERAGING) << "When returning the traction to the corner point, averaging of principal stresses is not supported\n";

    return CalculateCornerPoint();
}

Geo::PQ CoulombWithTensionCutOffImpl::ReturnStressAtCompressionCapZone(const Geo::PQ& rTrialPQ, const Matrix& rElasticMatrix) const
{
    const auto derivative_of_flow_function = mCompressionCapYieldSurface->DerivativeOfFlowFunction(rTrialPQ);
    const auto lambda = mCompressionCapYieldSurface->CalculatePlasticMultiplier(rTrialPQ, derivative_of_flow_function, rElasticMatrix);
    return Geo::PQ{rTrialPQ.Values() + lambda * prod(subrange(rElasticMatrix, 0, 2, 0, 2), derivative_of_flow_function)};
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
