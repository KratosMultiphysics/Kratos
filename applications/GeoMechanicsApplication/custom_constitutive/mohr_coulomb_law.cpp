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
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/coulomb_impl.h"
#include "custom_constitutive/mohr_coulomb_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

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
MohrCoulombLaw::MohrCoulombLaw()                                     = default;
MohrCoulombLaw::~MohrCoulombLaw()                                    = default;
MohrCoulombLaw::MohrCoulombLaw(MohrCoulombLaw&&) noexcept            = default;
MohrCoulombLaw& MohrCoulombLaw::operator=(MohrCoulombLaw&&) noexcept = default;

MohrCoulombLaw::MohrCoulombLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mpCoulombImpl(std::make_unique<CoulombImpl>())
{
}

ConstitutiveLaw::Pointer MohrCoulombLaw::Clone() const
{
    auto p_result           = std::make_shared<MohrCoulombLaw>(mpConstitutiveDimension->Clone());
    p_result->mStressVector = mStressVector;
    p_result->mStressVectorFinalized = mStressVectorFinalized;
    p_result->mStrainVectorFinalized = mStrainVectorFinalized;
    p_result->mpCoulombImpl          = mpCoulombImpl->Clone();
    p_result->mMaxRelativeOvershoot = mMaxRelativeOvershoot;
    p_result->mMaxNumberOfSubSteps =  mMaxNumberOfSubSteps;
    p_result->mCalculatedNumberOfSubSteps = mCalculatedNumberOfSubSteps;
    return p_result;
}

Vector& MohrCoulombLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

int& MohrCoulombLaw::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mpCoulombImpl->GetPlasticityStatus());
    }
    return rValue;
}

void MohrCoulombLaw::SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mStressVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType MohrCoulombLaw::WorkingSpaceDimension() { return mpConstitutiveDimension->GetDimension(); }

int MohrCoulombLaw::Check(const Properties&   rMaterialProperties,
                          const GeometryType& rElementGeometry,
                          const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    if (ConstitutiveLawUtilities::WantTensionCutOff(rMaterialProperties)) {
        check_properties.Check(
            GEO_TENSILE_STRENGTH,
            rMaterialProperties[GEO_COHESION] /
                std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    }
    check_properties.Check(YOUNG_MODULUS);
    constexpr auto max_value_poisson_ratio = 0.5;
    check_properties.Check(POISSON_RATIO, max_value_poisson_ratio);
    return result;
}

ConstitutiveLaw::StressMeasure MohrCoulombLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType MohrCoulombLaw::GetStrainSize() const { return mpConstitutiveDimension->GetStrainSize(); }

ConstitutiveLaw::StrainMeasure MohrCoulombLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool MohrCoulombLaw::IsIncremental() { return true; }

bool MohrCoulombLaw::RequiresInitializeMaterialResponse() { return true; }

void MohrCoulombLaw::InitializeMaterial(const Properties& rMaterialProperties, const Geometry<Node>&, const Vector&)
{
    mpCoulombImpl = std::make_unique<CoulombImpl>(rMaterialProperties);

    if (rMaterialProperties.Has(GEO_MAX_RELATIVE_OVERSHOOT)) {
        mMaxRelativeOvershoot = rMaterialProperties[GEO_MAX_RELATIVE_OVERSHOOT];
    }
    if (rMaterialProperties.Has(GEO_MAX_NUMBER_OF_SUB_STEPS)) {
        mMaxNumberOfSubSteps = rMaterialProperties[GEO_MAX_NUMBER_OF_SUB_STEPS];
    }
}

void MohrCoulombLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (!mIsModelInitialized) {
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }
}

void MohrCoulombLaw::GetLawFeatures(Features& rFeatures)
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

void MohrCoulombLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_properties = rParameters.GetMaterialProperties();

    if (rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rParameters.GetConstitutiveMatrix() =
            mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
    }
    if (!rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        return;
    }

    const auto elastic_matrix = mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);

    const Vector total_strain_increment = rParameters.GetStrainVector() - mStrainVectorFinalized;

    // Full elastic predictor over the entire strain increment.
    const Vector full_trial_stress_vector =
        mStressVectorFinalized + prod(elastic_matrix, total_strain_increment);
    const auto& [full_trial_principal_stresses, full_rotation_matrix] =
        StressStrainUtilities::CalculatePrincipalStressesAndRotationMatrix(full_trial_stress_vector);

    // If the whole step stays elastic, there is no need to sub-step.
    if (mpCoulombImpl->IsAdmissibleStressState(full_trial_principal_stresses)) {
        mStressVector                 = full_trial_stress_vector;
        rParameters.GetStressVector() = mStressVector;
        return;
    }

    std::size_t number_of_sub_steps = 1;
    if (mMaxNumberOfSubSteps > 1) {
        
        // only calculate the number of sub-steps once per solution step for stability
        if (mCalculatedNumberOfSubSteps == 0) {
            mCalculatedNumberOfSubSteps = CalculateAdaptiveNumberOfSubSteps(
            mpCoulombImpl, full_trial_principal_stresses, elastic_matrix);
        }
        
        number_of_sub_steps = mCalculatedNumberOfSubSteps;

    }

    // Running committed state for the sub-stepping (start from last finalized state)
    Vector committed_stress = mStressVectorFinalized;
    Vector committed_strain = mStrainVectorFinalized;

    for (std::size_t sub = 1; sub <= number_of_sub_steps; ++sub) {
        // strain at the end of this sub-step
        const Vector sub_strain =
            mStrainVectorFinalized +
            (static_cast<double>(sub) / static_cast<double>(number_of_sub_steps)) * total_strain_increment;

        // elastic predictor from the committed sub-step state
        const Vector trial_stress_vector =
            committed_stress + prod(elastic_matrix, sub_strain - committed_strain);

        const auto& [trial_principal_stresses, rotation_matrix] =
            StressStrainUtilities::CalculatePrincipalStressesAndRotationMatrix(trial_stress_vector);

        if (mpCoulombImpl->IsAdmissibleStressState(trial_principal_stresses)) {
            mStressVector = trial_stress_vector;
        } else {
            mpCoulombImpl->SaveKappaOfCoulombYieldSurface();
            auto mapped_principal_stresses = mpCoulombImpl->DoReturnMapping(
                trial_principal_stresses, elastic_matrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);

            if (const auto averaging_type = FindAveragingType(mapped_principal_stresses);
                averaging_type != Geo::PrincipalStresses::AveragingType::NO_AVERAGING) {
                const auto averaged_principal_trial_stress_vector =
                    AveragePrincipalStressComponents(trial_principal_stresses, averaging_type);
                mpCoulombImpl->RestoreKappaOfCoulombYieldSurface();
                mapped_principal_stresses = mpCoulombImpl->DoReturnMapping(
                    averaged_principal_trial_stress_vector, elastic_matrix, averaging_type);
                mapped_principal_stresses.Values()[1] =
                    mapped_principal_stresses.Values()[AveragingTypeToArrayIndex(averaging_type)];
            }
            mStressVector = StressStrainUtilities::RotatePrincipalStresses(
                mapped_principal_stresses.CopyTo<Vector>(), rotation_matrix,
                mpConstitutiveDimension->GetStrainSize());
        }

        committed_stress = mStressVector;
        committed_strain = sub_strain;
    }

    rParameters.GetStressVector() = mStressVector;
}

Vector MohrCoulombLaw::CalculateTrialStressVector(const Vector& rStrainVector, const Properties& rProperties) const
{
    return mStressVectorFinalized + prod(mpConstitutiveDimension->CalculateElasticConstitutiveTensor(rProperties),
                                         rStrainVector - mStrainVectorFinalized);
}

std::size_t MohrCoulombLaw::CalculateAdaptiveNumberOfSubSteps(std::unique_ptr<CoulombImpl>& rImpl,
                                                              const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                              const Matrix& rElasticMatrix)
{
    // make sure that kappa is not updated while calculating the number of required sub steps
    rImpl->SaveKappaOfCoulombYieldSurface();
    const auto mapped_principal_stresses = rImpl->DoReturnMapping(
        rTrialPrincipalStresses, rElasticMatrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
    rImpl->RestoreKappaOfCoulombYieldSurface();

    const Vector trial_values  = rTrialPrincipalStresses.CopyTo<Vector>();
    const Vector mapped_values = mapped_principal_stresses.CopyTo<Vector>();

    const auto overshoot          = norm_2(trial_values - mapped_values);
    const auto stress_scale       = std::max(norm_2(trial_values), 1.0e-12);
    const auto relative_overshoot = overshoot / stress_scale;

    const auto number_of_sub_steps =
        static_cast<int>(std::ceil(relative_overshoot / mMaxRelativeOvershoot));

    return std::clamp(number_of_sub_steps, 1, mMaxNumberOfSubSteps);
}

void MohrCoulombLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;

    mCalculatedNumberOfSubSteps = 0;
}

void MohrCoulombLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("mpCoulombImpl", mpCoulombImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
	rSerializer.save("MaxRelativeOvershoot", mMaxRelativeOvershoot);
	rSerializer.save("MaxNumberOfSubSteps", mMaxNumberOfSubSteps);
}

void MohrCoulombLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("mpCoulombImpl", mpCoulombImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
	rSerializer.load("MaxRelativeOvershoot", mMaxRelativeOvershoot);
	rSerializer.load("MaxNumberOfSubSteps", mMaxNumberOfSubSteps);
}

} // Namespace Kratos
