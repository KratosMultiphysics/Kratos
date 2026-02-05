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
//  Main authors:    Mohamed Nabi,
//                   Wijtze Pieter Kikstra
//

// Application includes
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>

namespace
{
using namespace Kratos;

std::size_t AveragingTypeToArrayIndex(CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    switch (AveragingType) {
        using enum CoulombYieldSurface::CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        return 0;
    case HIGHEST_PRINCIPAL_STRESSES:
        return 2;
    default:
        return 1;
    }
}

Vector AveragePrincipalStressComponents(const Vector& rPrincipalStressVector,
                                        CoulombYieldSurface::CoulombAveragingType AveragingType)
{
    auto result = rPrincipalStressVector;
    switch (AveragingType) {
        using enum CoulombYieldSurface::CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        std::fill(result.begin(), result.begin() + 1,
                  (rPrincipalStressVector[0] + rPrincipalStressVector[1]) * 0.5);
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        std::fill(result.begin() + 1, result.begin() + 2,
                  (rPrincipalStressVector[1] + rPrincipalStressVector[2]) * 0.5);
        break;
    default:
        break;
    }
    return result;
}

CoulombYieldSurface::CoulombAveragingType FindAveragingType(const Vector& rMappedPrincipalStressVector)
{
    using enum CoulombYieldSurface::CoulombAveragingType;
    if (rMappedPrincipalStressVector[0] < rMappedPrincipalStressVector[1]) {
        return LOWEST_PRINCIPAL_STRESSES;
    }
    if (rMappedPrincipalStressVector[1] < rMappedPrincipalStressVector[2]) {
        return HIGHEST_PRINCIPAL_STRESSES;
    }
    return NO_AVERAGING;
}

} // namespace

namespace Kratos
{

MohrCoulombWithTensionCutOff::MohrCoulombWithTensionCutOff(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

ConstitutiveLaw::Pointer MohrCoulombWithTensionCutOff::Clone() const
{
    auto p_result = std::make_shared<MohrCoulombWithTensionCutOff>(mpConstitutiveDimension->Clone());
    p_result->mStressVector                 = mStressVector;
    p_result->mStressVectorFinalized        = mStressVectorFinalized;
    p_result->mStrainVectorFinalized        = mStrainVectorFinalized;
    p_result->mCoulombWithTensionCutOffImpl = mCoulombWithTensionCutOffImpl;
    return p_result;
}

Vector& MohrCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

void MohrCoulombWithTensionCutOff::SetValue(const Variable<Vector>& rVariable,
                                            const Vector&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mStressVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType MohrCoulombWithTensionCutOff::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

int MohrCoulombWithTensionCutOff::Check(const Properties&   rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(
        GEO_TENSILE_STRENGTH,
        rMaterialProperties[GEO_COHESION] /
            std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    check_properties.Check(YOUNG_MODULUS);
    constexpr auto max_value_poisson_ratio = 0.5;
    check_properties.Check(POISSON_RATIO, max_value_poisson_ratio);
    return result;
}

ConstitutiveLaw::StressMeasure MohrCoulombWithTensionCutOff::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType MohrCoulombWithTensionCutOff::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

ConstitutiveLaw::StrainMeasure MohrCoulombWithTensionCutOff::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool MohrCoulombWithTensionCutOff::IsIncremental() { return true; }

bool MohrCoulombWithTensionCutOff::RequiresInitializeMaterialResponse() { return true; }

void MohrCoulombWithTensionCutOff::InitializeMaterial(const Properties& rMaterialProperties,
                                                      const Geometry<Node>&,
                                                      const Vector&)
{
    mCoulombWithTensionCutOffImpl = CoulombWithTensionCutOffImpl{rMaterialProperties};
}

void MohrCoulombWithTensionCutOff::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (!mIsModelInitialized) {
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }
}

void MohrCoulombWithTensionCutOff::GetLawFeatures(Features& rFeatures)
{
    auto options = Flags{};
    options.Set(mpConstitutiveDimension->GetSpatialType());
    options.Set(ConstitutiveLaw::INFINITESIMAL_STRAINS);
    options.Set(ConstitutiveLaw::ISOTROPIC);
    rFeatures.SetOptions(options);

    rFeatures.SetStrainMeasure(StrainMeasure_Infinitesimal);
    rFeatures.SetStrainSize(GetStrainSize());
    rFeatures.SetSpaceDimension(WorkingSpaceDimension());
}

void MohrCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_properties = rParameters.GetMaterialProperties();

    if (rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rParameters.GetConstitutiveMatrix() = mpConstitutiveDimension->CalculateElasticMatrix(r_properties);
    }

    const auto trial_stress_vector = CalculateTrialStressVector(rParameters.GetStrainVector(), r_properties);

    Vector principal_trial_stress_vector;
    Matrix rotation_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(
        trial_stress_vector, principal_trial_stress_vector, rotation_matrix);

    if (auto trial_sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(principal_trial_stress_vector);
        mCoulombWithTensionCutOffImpl.IsAdmissibleSigmaTau(trial_sigma_tau)) {
        mStressVector = trial_stress_vector;
    } else {
        mCoulombWithTensionCutOffImpl.SaveKappaOfCoulombYieldSurface();
        auto mapped_sigma_tau = mCoulombWithTensionCutOffImpl.DoReturnMapping(
            trial_sigma_tau, CoulombYieldSurface::CoulombAveragingType::NO_AVERAGING);
        auto mapped_principal_stress_vector = StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
            mapped_sigma_tau, principal_trial_stress_vector);

        // for interchanging principal stresses, retry mapping with averaged principal stresses.

        if (const auto averaging_type = FindAveragingType(mapped_principal_stress_vector);
            averaging_type != CoulombYieldSurface::CoulombAveragingType::NO_AVERAGING) {
            const auto averaged_principal_trial_stress_vector =
                AveragePrincipalStressComponents(principal_trial_stress_vector, averaging_type);
            trial_sigma_tau = StressStrainUtilities::TransformPrincipalStressesToSigmaTau(
                averaged_principal_trial_stress_vector);
            mCoulombWithTensionCutOffImpl.RestoreKappaOfCoulombYieldSurface();
            mapped_sigma_tau = mCoulombWithTensionCutOffImpl.DoReturnMapping(trial_sigma_tau, averaging_type);
            mapped_principal_stress_vector = StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
                mapped_sigma_tau, averaged_principal_trial_stress_vector);
            mapped_principal_stress_vector[1] =
                mapped_principal_stress_vector[AveragingTypeToArrayIndex(averaging_type)];
        }
        mStressVector = StressStrainUtilities::RotatePrincipalStresses(
            mapped_principal_stress_vector, rotation_matrix, mpConstitutiveDimension->GetStrainSize());
    }

    rParameters.GetStressVector() = mStressVector;
}

Vector MohrCoulombWithTensionCutOff::CalculateTrialStressVector(const Vector&     rStrainVector,
                                                                const Properties& rProperties) const
{
    return mStressVectorFinalized + prod(mpConstitutiveDimension->CalculateElasticMatrix(rProperties),
                                         rStrainVector - mStrainVectorFinalized);
}

void MohrCoulombWithTensionCutOff::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void MohrCoulombWithTensionCutOff::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void MohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos
