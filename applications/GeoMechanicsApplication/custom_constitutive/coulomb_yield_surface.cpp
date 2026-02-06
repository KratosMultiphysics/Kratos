// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/geo_sigma_tau.hpp"
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/function_object_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/string_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <cmath>

namespace
{

using namespace Kratos;

std::string GetCoulombHardeningTypeFrom(const Properties& rMaterialProperties)
{
    return GeoStringUtilities::ToLower(rMaterialProperties[GEO_COULOMB_HARDENING_TYPE]);
}

Geo::KappaDependentFunction MakeFrictionAngleCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return FunctionObjectUtilities::MakeConstantFunction(
            ConstitutiveLawUtilities::GetFrictionAngleInRadians(rMaterialProperties));
    }

    if (hardening_type == "linear") {
        return FunctionObjectUtilities::MakeLinearFunction(
            ConstitutiveLawUtilities::GetFrictionAngleInRadians(rMaterialProperties),
            rMaterialProperties[GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS][0]);
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the friction angle of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

Geo::KappaDependentFunction MakeCohesionCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return FunctionObjectUtilities::MakeConstantFunction(
            ConstitutiveLawUtilities::GetCohesion(rMaterialProperties));
    }

    if (hardening_type == "linear") {
        return FunctionObjectUtilities::MakeLinearFunction(
            ConstitutiveLawUtilities::GetCohesion(rMaterialProperties),
            rMaterialProperties[GEO_COHESION_FUNCTION_COEFFICIENTS][0]);
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cohesion of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

Geo::KappaDependentFunction MakeDilatancyAngleCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return FunctionObjectUtilities::MakeConstantFunction(
            MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]));
    }

    if (hardening_type == "linear") {
        return FunctionObjectUtilities::MakeLinearFunction(
            MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]),
            rMaterialProperties[GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS][0]);
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the dilatancy angle of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

} // namespace

namespace Kratos
{

CoulombYieldSurface::CoulombYieldSurface()
{
    mMaterialProperties.SetValue(GEO_COULOMB_HARDENING_TYPE, "None");
    mMaterialProperties.SetValue(GEO_FRICTION_ANGLE, 0.0);
    mMaterialProperties.SetValue(GEO_COHESION, 0.0);
    mMaterialProperties.SetValue(GEO_DILATANCY_ANGLE, 0.0);

    InitializeKappaDependentFunctions();
}

CoulombYieldSurface::CoulombYieldSurface(const Properties& rMaterialProperties)
    : mMaterialProperties{rMaterialProperties}
{
    // For backward compatibility, if no hardening type is given, we assume no hardening at all
    if (!mMaterialProperties.Has(GEO_COULOMB_HARDENING_TYPE)) {
        mMaterialProperties.SetValue(GEO_COULOMB_HARDENING_TYPE, "None");
    }

    CheckMaterialProperties();
    InitializeKappaDependentFunctions();
}

double CoulombYieldSurface::GetFrictionAngleInRadians() const
{
    return mFrictionAngleCalculator(mKappa);
}

double CoulombYieldSurface::GetCohesion() const { return mCohesionCalculator(mKappa); }

double CoulombYieldSurface::GetDilatancyAngleInRadians() const
{
    return mDilatancyAngleCalculator(mKappa);
}

double CoulombYieldSurface::GetKappa() const { return mKappa; }

void CoulombYieldSurface::SetKappa(double kappa) { mKappa = kappa; }

// At some point in time we would like to get rid of this API. For now, just forward the request.
double CoulombYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return YieldFunctionValue(Geo::SigmaTau{rSigmaTau});
}

double CoulombYieldSurface::YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const
{
    return rSigmaTau.Tau() + rSigmaTau.Sigma() * std::sin(GetFrictionAngleInRadians()) -
           GetCohesion() * std::cos(GetFrictionAngleInRadians());
}

double CoulombYieldSurface::YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return YieldFunctionValue(StressStrainUtilities::TransformPrincipalStressesToSigmaTau(rPrincipalStresses));
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return DerivativeOfFlowFunction(Geo::SigmaTau{rSigmaTau}, CoulombAveragingType::NO_AVERAGING);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Geo::SigmaTau&, CoulombAveragingType AveragingType) const
{
    const auto sin_psi = std::sin(GetDilatancyAngleInRadians());
    switch (AveragingType) {
        using enum CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        return UblasUtilities::CreateVector({-(1.0 - 3.0 * sin_psi) / 4.0, (3.0 - sin_psi) / 4.0});
    case NO_AVERAGING:
        return UblasUtilities::CreateVector({sin_psi, 1.0});
    case HIGHEST_PRINCIPAL_STRESSES:
        return UblasUtilities::CreateVector({(1.0 + 3.0 * sin_psi) / 4.0, (3.0 + sin_psi) / 4.0});
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
}

double CoulombYieldSurface::CalculateApex() const
{
    return GetCohesion() / std::tan(GetFrictionAngleInRadians());
}

void CoulombYieldSurface::InitializeKappaDependentFunctions()
{
    mFrictionAngleCalculator  = MakeFrictionAngleCalculator(mMaterialProperties);
    mCohesionCalculator       = MakeCohesionCalculator(mMaterialProperties);
    mDilatancyAngleCalculator = MakeDilatancyAngleCalculator(mMaterialProperties);
}

double CoulombYieldSurface::CalculatePlasticMultiplier(const Geo::SigmaTau& rSigmaTau,
                                                       const Vector& rDerivativeOfFlowFunction) const
{
    const auto sin_phi   = std::sin(GetFrictionAngleInRadians());
    const auto numerator = sin_phi * rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1];
    return (GetCohesion() * std::cos(GetFrictionAngleInRadians()) - rSigmaTau.Sigma() * sin_phi -
            rSigmaTau.Tau()) /
           numerator;
}

double CoulombYieldSurface::CalculateEquivalentPlasticStrainIncrement(const Geo::SigmaTau& rSigmaTau,
                                                                      CoulombAveragingType AveragingType) const
{
    const auto derivative              = DerivativeOfFlowFunction(rSigmaTau, AveragingType);
    const auto principal_stress_vector = UblasUtilities::CreateVector(
        {(derivative[0] + derivative[1]) / 2.0, 0.0, (derivative[0] - derivative[1]) / 2.0});
    const auto mean = std::accumulate(principal_stress_vector.begin(), principal_stress_vector.end(), 0.0) /
                      static_cast<double>(principal_stress_vector.size());
    auto deviatoric_principle_stress_vector = Vector{3};
    std::ranges::transform(principal_stress_vector, deviatoric_principle_stress_vector.begin(),
                           [mean](auto sigma) { return sigma - mean; });
    return -std::sqrt(2.0 / 3.0) * MathUtils<>::Norm(deviatoric_principle_stress_vector) *
           CalculatePlasticMultiplier(rSigmaTau, DerivativeOfFlowFunction(rSigmaTau, AveragingType));
}

void CoulombYieldSurface::CheckMaterialProperties() const
{
    const CheckProperties check_properties(mMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COHESION);
    check_properties.Check(GEO_FRICTION_ANGLE);
    check_properties.Check(GEO_DILATANCY_ANGLE, mMaterialProperties[GEO_FRICTION_ANGLE]);

    if (GetCoulombHardeningTypeFrom(mMaterialProperties) != "none") {
        CheckHardeningCoefficients(GEO_COHESION_FUNCTION_COEFFICIENTS, check_properties);
        CheckHardeningCoefficients(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, check_properties);
        CheckHardeningCoefficients(GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS, check_properties);
    }
}

void CoulombYieldSurface::CheckHardeningCoefficients(const Variable<Vector>& rCoefficientsVariable,
                                                     const CheckProperties&  rChecker) const
{
    rChecker.CheckAvailability(rCoefficientsVariable);

    const auto& coefficients = mMaterialProperties[rCoefficientsVariable];
    for (auto i = std::size_t{0}; i < coefficients.size(); ++i) {
        KRATOS_ERROR_IF(coefficients[i] < 0.0) << "Entry " << i << " in " << rCoefficientsVariable.Name()
                                               << " out of range. Value: " << coefficients[i];
    }
}

void CoulombYieldSurface::save(Serializer& rSerializer) const
{
    rSerializer.save("Kappa", mKappa);
    rSerializer.save("MaterialProperties", mMaterialProperties);

    // Members `mFrictionAngleCalculator`, `mCohesionCalculator`, and `mDilatancyAngleCalculator`
    // will be reconstructed using data from `mMaterialProperties`. No additional serialization is
    // needed.
}

void CoulombYieldSurface::load(Serializer& rSerializer)
{
    rSerializer.load("Kappa", mKappa);
    rSerializer.load("MaterialProperties", mMaterialProperties);

    InitializeKappaDependentFunctions();
}

} // Namespace Kratos
