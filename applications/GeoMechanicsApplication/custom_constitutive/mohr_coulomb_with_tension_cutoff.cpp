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
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>

namespace
{

using namespace Kratos;

void CheckProperty(const Properties&       rMaterialProperties,
                   const Variable<double>& rVariable,
                   std::optional<double>   MaxValue = std::nullopt)
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(rVariable))
        << rVariable.Name() << " is not defined for property " << rMaterialProperties.Id() << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[rVariable] < 0.0 ||
                    (MaxValue.has_value() && rMaterialProperties[rVariable] > MaxValue.value()))
        << "value of " << rVariable.Name() << " for property " << rMaterialProperties.Id()
        << " is out of range: " << rMaterialProperties[rVariable] << " is not in [0.0, "
        << (MaxValue ? std::to_string(*MaxValue) + "]" : "->") << std::endl;
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
    p_result->mStressVector          = mStressVector;
    p_result->mStressVectorFinalized = mStressVectorFinalized;
    p_result->mStrainVectorFinalized = mStrainVectorFinalized;
    p_result->mCoulombYieldSurface   = mCoulombYieldSurface;
    p_result->mTensionCutOff         = mTensionCutOff;
    return p_result;
}

Vector& MohrCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rThisVariable, rValue);
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

    CheckProperty(rMaterialProperties, GEO_COHESION);
    CheckProperty(rMaterialProperties, GEO_FRICTION_ANGLE);
    CheckProperty(rMaterialProperties, GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
    CheckProperty(rMaterialProperties, GEO_TENSILE_STRENGTH,
                  rMaterialProperties[GEO_COHESION] /
                      std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    CheckProperty(rMaterialProperties, YOUNG_MODULUS);
    CheckProperty(rMaterialProperties, POISSON_RATIO, 0.5);
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
    mCoulombYieldSurface =
        CoulombYieldSurface(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]),
                            rMaterialProperties[GEO_COHESION],
                            MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]));
    mTensionCutOff = TensionCutoff(rMaterialProperties[GEO_TENSILE_STRENGTH]);
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
    const auto& r_prop = rParameters.GetMaterialProperties();

    if (rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rParameters.GetConstitutiveMatrix() =
            mpConstitutiveDimension->CalculateElasticMatrix(r_prop[YOUNG_MODULUS], r_prop[POISSON_RATIO]);
    }

    const auto trail_stress_vector = CalculateTrialStressVector(
        rParameters.GetStrainVector(), r_prop[YOUNG_MODULUS], r_prop[POISSON_RATIO]);

    Vector principal_trial_stress_vector;
    Matrix rotation_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(
        trail_stress_vector, principal_trial_stress_vector, rotation_matrix);

    auto trial_sigma_tau =
        StressStrainUtilities::TransformPrincipalStressesToSigmaAndTau(principal_trial_stress_vector);

    while (!ConstitutiveLawUtilities::IsAdmissibleSigmaTauStressState(
        trial_sigma_tau, mCoulombYieldSurface, mTensionCutOff)) {
        trial_sigma_tau = ConstitutiveLawUtilities::ReturnMappingToCoulombWithTensionCutOff(
            r_prop, trial_sigma_tau, mCoulombYieldSurface, mTensionCutOff);
        principal_trial_stress_vector = StressStrainUtilities::TransformSigmaAndTauToPrincipalStresses(
            trial_sigma_tau, principal_trial_stress_vector);

        StressStrainUtilities::ReorderEigenValuesAndVectors(principal_trial_stress_vector, rotation_matrix);
        trial_sigma_tau =
            StressStrainUtilities::TransformPrincipalStressesToSigmaAndTau(principal_trial_stress_vector);
    }

    mStressVector = StressStrainUtilities::RotatePrincipalStresses(
        principal_trial_stress_vector, rotation_matrix, mpConstitutiveDimension->GetStrainSize());

    rParameters.GetStressVector() = mStressVector;
}

Vector MohrCoulombWithTensionCutOff::CalculateTrialStressVector(const Vector& rStrainVector,
                                                                double        YoungsModulus,
                                                                double        PoissonsRatio) const
{
    return mStressVectorFinalized +
           prod(mpConstitutiveDimension->CalculateElasticMatrix(YoungsModulus, PoissonsRatio),
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
    rSerializer.save("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.save("TensionCutOff", mTensionCutOff);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void MohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.load("TensionCutOff", mTensionCutOff);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}
} // Namespace Kratos
