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
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <algorithm>
#include <cmath>
#include <limits>

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

int& MohrCoulombWithTensionCutOff::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mCoulombWithTensionCutOffImpl.GetPlasticityStatus());
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

    const auto calculate_tangent = rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const auto calculate_stress  = rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS);

    if (calculate_tangent && calculate_stress) {
        rParameters.GetConstitutiveMatrix() =
            CalculateNumericalTangentMatrix(rParameters.GetStrainVector(), r_properties);
    } else if (calculate_tangent) {
        rParameters.GetConstitutiveMatrix() =
            mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
    }
    if (!calculate_stress) {
        return;
    }

    mStressVector =
        CalculateStressVector(rParameters.GetStrainVector(), r_properties, mCoulombWithTensionCutOffImpl);
    rParameters.GetStressVector() = mStressVector;
}

Vector MohrCoulombWithTensionCutOff::CalculateTrialStressVector(const Vector&     rStrainVector,
                                                                const Properties& rProperties) const
{
    return mStressVectorFinalized + prod(mpConstitutiveDimension->CalculateElasticConstitutiveTensor(rProperties),
                                         rStrainVector - mStrainVectorFinalized);
}

Vector MohrCoulombWithTensionCutOff::CalculateStressVector(
    const Vector&                rStrainVector,
    const Properties&            rProperties,
    CoulombWithTensionCutOffImpl& rCoulombWithTensionCutOffImpl) const
{
    const auto elastic_constitutive_tensor = mpConstitutiveDimension->CalculateElasticConstitutiveTensor(rProperties);
    const auto trial_stress_vector         = CalculateTrialStressVector(rStrainVector, rProperties);
    const auto& [trial_principal_stresses, rotation_matrix] =
        StressStrainUtilities::CalculatePrincipalStressesAndRotationMatrix(trial_stress_vector);

    if (rCoulombWithTensionCutOffImpl.IsAdmissibleStressState(trial_principal_stresses)) {
        return trial_stress_vector;
    } else {
        rCoulombWithTensionCutOffImpl.SaveKappaOfCoulombYieldSurface();
        auto mapped_principal_stresses = rCoulombWithTensionCutOffImpl.DoReturnMapping(
            trial_principal_stresses, elastic_constitutive_tensor,
            Geo::PrincipalStresses::AveragingType::NO_AVERAGING);

        // For interchanging principal stresses, retry mapping with averaged principal stresses.
        if (const auto averaging_type = FindAveragingType(mapped_principal_stresses);
            averaging_type != Geo::PrincipalStresses::AveragingType::NO_AVERAGING) {
            const auto averaged_principal_trial_stress_vector =
                AveragePrincipalStressComponents(trial_principal_stresses, averaging_type);
            rCoulombWithTensionCutOffImpl.RestoreKappaOfCoulombYieldSurface();
            mapped_principal_stresses = rCoulombWithTensionCutOffImpl.DoReturnMapping(
                averaged_principal_trial_stress_vector, elastic_constitutive_tensor, averaging_type);
            mapped_principal_stresses.Values()[1] =
                mapped_principal_stresses.Values()[AveragingTypeToArrayIndex(averaging_type)];
        }
        return StressStrainUtilities::RotatePrincipalStresses(
            mapped_principal_stresses.CopyTo<Vector>(), rotation_matrix,
            mpConstitutiveDimension->GetStrainSize());
    }
}

Matrix MohrCoulombWithTensionCutOff::CalculateNumericalTangentMatrix(const Vector&     rStrainVector,
                                                                     const Properties& rProperties) const
{
    constexpr auto perturbation_coefficient_1 = 1.0e-5;
    constexpr auto perturbation_coefficient_2 = 1.0e-10;
    constexpr auto perturbation_threshold     = 1.0e-8;

    const auto strain_size = rStrainVector.size();
    Matrix     result      = ZeroMatrix(strain_size, strain_size);
    const auto minmax_strain =
        std::minmax_element(rStrainVector.begin(), rStrainVector.end(), [](double Left, double Right) {
            return std::abs(Left) < std::abs(Right);
        });
    const auto min_abs_strain = std::abs(*minmax_strain.first);
    const auto max_abs_strain = std::abs(*minmax_strain.second);

    for (auto i = SizeType{0}; i < strain_size; ++i) {
        const auto strain_based_perturbation =
            perturbation_coefficient_1 *
            (std::abs(rStrainVector[i]) > std::numeric_limits<double>::epsilon()
                 ? std::abs(rStrainVector[i])
                 : min_abs_strain);
        const auto perturbation = std::max({perturbation_threshold, strain_based_perturbation,
                                            perturbation_coefficient_2 * max_abs_strain});

        auto strain_plus  = rStrainVector;
        auto strain_minus = rStrainVector;
        strain_plus[i] += perturbation;
        strain_minus[i] -= perturbation;

        auto coulomb_with_tension_cut_off_impl_plus  = mCoulombWithTensionCutOffImpl;
        auto coulomb_with_tension_cut_off_impl_minus = mCoulombWithTensionCutOffImpl;
        const auto stress_plus =
            CalculateStressVector(strain_plus, rProperties, coulomb_with_tension_cut_off_impl_plus);
        const auto stress_minus =
            CalculateStressVector(strain_minus, rProperties, coulomb_with_tension_cut_off_impl_minus);

        column(result, i) = (stress_plus - stress_minus) / (2.0 * perturbation);
    }

    return result;
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
