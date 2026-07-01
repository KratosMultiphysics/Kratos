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
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff_elastoplastic_tangent_matrix.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "custom_utilities/constitutive_tangent_operator_calculator.h"

#include <cmath>

namespace
{
using namespace Kratos;

std::size_t AveragingTypeToArrayIndex(Geo::PrincipalStresses::AveragingType AveragingType)
{
    switch (AveragingType) {
        using enum Geo::PrincipalStresses::AveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        return 0;
    case HIGHEST_PRINCIPAL_STRESSES:
        return 2;
    default:
        return 1;
    }
}

Geo::PrincipalStresses AveragePrincipalStressComponents(const Geo::PrincipalStresses& rPrincipalStressVector,
                                                        Geo::PrincipalStresses::AveragingType AveragingType)
{
    auto result = rPrincipalStressVector;
    switch (AveragingType) {
        using enum Geo::PrincipalStresses::AveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        std::fill(result.Values().begin(), result.Values().begin() + std::ptrdiff_t{2},
                  (rPrincipalStressVector.Values()[0] + rPrincipalStressVector.Values()[1]) * 0.5);
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        std::fill(result.Values().end() - std::ptrdiff_t{2}, result.Values().end(),
                  (rPrincipalStressVector.Values()[1] + rPrincipalStressVector.Values()[2]) * 0.5);
        break;
    default:
        break;
    }
    return result;
}

Geo::PrincipalStresses::AveragingType FindAveragingType(const Geo::PrincipalStresses& rMappedPrincipalStressVector)
{
    using enum Geo::PrincipalStresses::AveragingType;
    if (rMappedPrincipalStressVector.Values()[0] < rMappedPrincipalStressVector.Values()[1]) {
        return LOWEST_PRINCIPAL_STRESSES;
    }
    if (rMappedPrincipalStressVector.Values()[1] < rMappedPrincipalStressVector.Values()[2]) {
        return HIGHEST_PRINCIPAL_STRESSES;
    }
    return NO_AVERAGING;
}

} // namespace

namespace Kratos
{

MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

ConstitutiveLaw::Pointer MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::Clone() const
{
    auto p_result = std::make_shared<MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix>(mpConstitutiveDimension->Clone());
    p_result->mStressVector                 = mStressVector;
    p_result->mStressVectorFinalized        = mStressVectorFinalized;
    p_result->mStrainVectorFinalized        = mStrainVectorFinalized;
    p_result->mCoulombWithTensionCutOffImpl = mCoulombWithTensionCutOffImpl;
    p_result->mIsModelInitialized           = mIsModelInitialized;
    return p_result;
}

Vector& MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

int& MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mCoulombWithTensionCutOffImpl.GetPlasticityStatus());
    }
    return rValue;
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::SetValue(const Variable<Vector>& rVariable,
                                            const Vector&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mStressVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

int MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::Check(const Properties&   rMaterialProperties,
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

ConstitutiveLaw::StressMeasure MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

ConstitutiveLaw::StrainMeasure MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::IsIncremental() { return true; }

bool MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::RequiresInitializeMaterialResponse() { return true; }

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::InitializeMaterial(const Properties& rMaterialProperties,
                                                      const Geometry<Node>&,
                                                      const Vector&)
{
    mCoulombWithTensionCutOffImpl = CoulombImpl{rMaterialProperties};
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (!mIsModelInitialized) {
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::GetLawFeatures(Features& rFeatures)
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

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_properties = rParameters.GetMaterialProperties();
    auto& r_options = rParameters.GetOptions();
    const bool compute_stress = r_options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool compute_constitutive_tensor =
        r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool use_numerical_tangent_operator =
        r_properties.Has(GEO_USE_NUMERICAL_TANGENT_OPERATOR)
            ? r_properties[GEO_USE_NUMERICAL_TANGENT_OPERATOR]
            : false;

    if (compute_constitutive_tensor && !use_numerical_tangent_operator) {
        rParameters.GetConstitutiveMatrix() =
            mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
    }

    ConstitutiveLaw::Pointer p_tangent_base_law;
    if (compute_constitutive_tensor && use_numerical_tangent_operator) {
        p_tangent_base_law = Clone();
    }

    if (!compute_stress) {
        if (compute_constitutive_tensor && use_numerical_tangent_operator) {
            const Vector unperturbed_stress = rParameters.GetStressVector();
            auto p_unperturbed_response_law = p_tangent_base_law->Clone();

            r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
            r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            p_unperturbed_response_law->CalculateMaterialResponseCauchy(rParameters);
            r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            ConstitutiveTangentOperatorCalculator::CalculateTangentTensor(
                rParameters, *p_tangent_base_law);
            if (norm_frobenius(rParameters.GetConstitutiveMatrix()) < 1.0e-8) {
                constexpr double regularization_factor = 1.0;
                rParameters.GetConstitutiveMatrix() =
                    regularization_factor *
                    mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
            }
            rParameters.GetStressVector() = unperturbed_stress;
        }
        return;
    }

    const auto trial_stress_vector = CalculateTrialStressVector(rParameters.GetStrainVector(), r_properties);
    const auto& [trial_principal_stresses, rotation_matrix] =
        StressStrainUtilities::CalculatePrincipalStressesAndRotationMatrix(trial_stress_vector);

    if (mCoulombWithTensionCutOffImpl.IsAdmissibleStressState(trial_principal_stresses)) {
        mStressVector = trial_stress_vector;
    } else {
        mCoulombWithTensionCutOffImpl.SaveKappaOfCoulombYieldSurface();
        auto mapped_principal_stresses = mCoulombWithTensionCutOffImpl.DoReturnMapping(
            trial_principal_stresses, mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties),
            Geo::PrincipalStresses::AveragingType::NO_AVERAGING);

        // For interchanging principal stresses, retry mapping with averaged principal stresses.
        if (const auto averaging_type = FindAveragingType(mapped_principal_stresses);
            averaging_type != Geo::PrincipalStresses::AveragingType::NO_AVERAGING) {
            const auto averaged_principal_trial_stress_vector =
                AveragePrincipalStressComponents(trial_principal_stresses, averaging_type);
            mCoulombWithTensionCutOffImpl.RestoreKappaOfCoulombYieldSurface();
            mapped_principal_stresses = mCoulombWithTensionCutOffImpl.DoReturnMapping(
                averaged_principal_trial_stress_vector,
                mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties), averaging_type);
            mapped_principal_stresses.Values()[1] =
                mapped_principal_stresses.Values()[AveragingTypeToArrayIndex(averaging_type)];
        }
        mStressVector = StressStrainUtilities::RotatePrincipalStresses(
            mapped_principal_stresses.CopyTo<Vector>(), rotation_matrix,
            mpConstitutiveDimension->GetStrainSize());
    }
    rParameters.GetStressVector() = mStressVector;

    if (compute_constitutive_tensor && use_numerical_tangent_operator) {
        ConstitutiveTangentOperatorCalculator::CalculateTangentTensor(
            rParameters, *p_tangent_base_law);
        if (norm_frobenius(rParameters.GetConstitutiveMatrix()) < 1.0e-8) {
            constexpr double regularization_factor = 1.0;
            rParameters.GetConstitutiveMatrix() =
                regularization_factor *
                mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
        }
    }
}

Vector MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::CalculateTrialStressVector(const Vector&     rStrainVector,
                                                                const Properties& rProperties) const
{
    return mStressVectorFinalized + prod(mpConstitutiveDimension->CalculateElasticConstitutiveTensor(rProperties),
                                         rStrainVector - mStrainVectorFinalized);
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void MohrCoulombWithTensionCutOffElastoPlasticTangentMatrix::load(Serializer& rSerializer)
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
