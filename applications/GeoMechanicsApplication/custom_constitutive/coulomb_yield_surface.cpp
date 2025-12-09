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
//

#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace
{

using namespace Kratos;

CoulombYieldSurface::KappaDependentFunction MakeConstantFunction(double Value)
{
    return [Value](double /* unused kappa */) { return Value; };
}

CoulombYieldSurface::KappaDependentFunction MakeLinearFunction(double Value, double Coefficient)
{
    return [Value, Coefficient](double Kappa) { return Value + Coefficient * Kappa; };
}

std::string GetCoulombHardeningTypeFrom(const Properties& rMaterialProperties)
{
    auto result   = rMaterialProperties[GEO_COULOMB_HARDENING_TYPE];
    auto to_lower = [](auto character) { return std::tolower(character); };
    std::ranges::transform(result, result.begin(), to_lower);
    return result;
}

CoulombYieldSurface::KappaDependentFunction MakeFrictionAngleCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantFunction(ConstitutiveLawUtilities::GetFrictionAngleInRadians(rMaterialProperties));
    }

    if (hardening_type == "linear") {
        return MakeLinearFunction(ConstitutiveLawUtilities::GetFrictionAngleInRadians(rMaterialProperties),
                                  rMaterialProperties[GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS][0]);
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the friction angle of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

CoulombYieldSurface::KappaDependentFunction MakeCohesionCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantFunction(ConstitutiveLawUtilities::GetCohesion(rMaterialProperties));
    }

    if (hardening_type == "linear") {
        return MakeLinearFunction(ConstitutiveLawUtilities::GetCohesion(rMaterialProperties),
                                  rMaterialProperties[GEO_COHESION_FUNCTION_COEFFICIENTS][0]);
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cohesion of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

CoulombYieldSurface::KappaDependentFunction MakeDilatancyAngleCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCoulombHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantFunction(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]));
    }

    if (hardening_type == "linear") {
        return MakeLinearFunction(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]),
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

double CoulombYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return rSigmaTau[1] + rSigmaTau[0] * std::sin(GetFrictionAngleInRadians()) -
           GetCohesion() * std::cos(GetFrictionAngleInRadians());
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return DerivativeOfFlowFunction(rSigmaTau, CoulombAveragingType::NO_AVERAGING);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const
{
    const auto sin_psi = std::sin(GetDilatancyAngleInRadians());
    Vector     result(2);
    switch (AveragingType) {
        using enum CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        result <<= -(1.0 - 3.0 * sin_psi) / 4.0, (3.0 - sin_psi) / 4.0;
        break;
    case NO_AVERAGING:
        result <<= sin_psi, 1.0;
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        result <<= (1.0 + 3.0 * sin_psi) / 4.0, (3.0 + sin_psi) / 4.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
    return result;
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

double CoulombYieldSurface::CalculatePlasticMultiplier(const Vector& rSigmaTau,
                                                       const Vector& rDerivativeOfFlowFunction) const
{
    const auto sin_phi   = std::sin(GetFrictionAngleInRadians());
    const auto numerator = sin_phi * rDerivativeOfFlowFunction[0] + rDerivativeOfFlowFunction[1];
    return (GetCohesion() * std::cos(GetFrictionAngleInRadians()) - rSigmaTau[0] * sin_phi - rSigmaTau[1]) / numerator;
}

double CoulombYieldSurface::CalculateEquivalentPlasticStrainIncrement(const Vector& rSigmaTau,
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
