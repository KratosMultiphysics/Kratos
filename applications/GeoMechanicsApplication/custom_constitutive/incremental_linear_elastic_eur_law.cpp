// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_constitutive/incremental_linear_elastic_eur_law.h"

#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/mat_variables.h"
#include "utilities/math_utils.h"

namespace
{

using namespace Kratos;

void SetEntriesAboveDiagonalToZero(Matrix& rMatrix)
{
    for (auto i = std::size_t{0}; i < rMatrix.size1() - 1; ++i) {
        for (auto j = i + 1; j < rMatrix.size2(); ++j) {
            rMatrix(i, j) = 0.0;
        }
    }
}

void SetEntriesBelowDiagonalToZero(Matrix& rMatrix)
{
    for (auto i = std::size_t{1}; i < rMatrix.size1(); ++i) {
        for (auto j = std::size_t{0}; j < i; ++j) {
            rMatrix(i, j) = 0.0;
        }
    }
}

void SetShearEntriesToZero(Matrix& rMatrix, std::size_t NumberOfNormalComponents)
{
    for (auto i = NumberOfNormalComponents; i < rMatrix.size1(); ++i) {
        rMatrix(i, i) = 0.0;
    }
}

} // namespace

namespace Kratos
{

GeoIncrementalLinearElasticEurLaw::GeoIncrementalLinearElasticEurLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : GeoLinearElasticLaw{},
      mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

GeoIncrementalLinearElasticEurLaw::GeoIncrementalLinearElasticEurLaw(const GeoIncrementalLinearElasticEurLaw& rOther)
    : GeoLinearElasticLaw(rOther),
      mStressVector(rOther.mStressVector),
      mStressVectorFinalized(rOther.mStressVectorFinalized),
      mDeltaStrainVector(rOther.mDeltaStrainVector),
      mStrainVectorFinalized(rOther.mStrainVectorFinalized),
      mIsModelInitialized(rOther.mIsModelInitialized)
{
    if (rOther.mpConstitutiveDimension) {
        mpConstitutiveDimension = rOther.mpConstitutiveDimension->Clone();
    }
}

GeoIncrementalLinearElasticEurLaw& GeoIncrementalLinearElasticEurLaw::operator=(const GeoIncrementalLinearElasticEurLaw& rOther)
{
    GeoLinearElasticLaw::operator=(rOther);
    mStressVector          = rOther.mStressVector;
    mStressVectorFinalized = rOther.mStressVectorFinalized;
    mDeltaStrainVector     = rOther.mDeltaStrainVector;
    mStrainVectorFinalized = rOther.mStrainVectorFinalized;
    mIsModelInitialized    = rOther.mIsModelInitialized;

    if (rOther.mpConstitutiveDimension) {
        mpConstitutiveDimension = rOther.mpConstitutiveDimension->Clone();
    }

    return *this;
}

ConstitutiveLaw::Pointer GeoIncrementalLinearElasticEurLaw::Clone() const
{
    return Kratos::make_shared<GeoIncrementalLinearElasticEurLaw>(*this);
}

bool& GeoIncrementalLinearElasticEurLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }
    return rValue;
}

void GeoIncrementalLinearElasticEurLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(mpConstitutiveDimension->GetSpatialType());
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize     = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

GeoIncrementalLinearElasticEurLaw::SizeType GeoIncrementalLinearElasticEurLaw::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

GeoIncrementalLinearElasticEurLaw::SizeType GeoIncrementalLinearElasticEurLaw::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

bool GeoIncrementalLinearElasticEurLaw::IsIncremental() { return true; }

int GeoIncrementalLinearElasticEurLaw::Check(const Properties&   rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const ProcessInfo&  rCurrentProcessInfo) const
{
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "parameters of material",
                                           CheckProperties::Bounds::AllExclusive);
    check_properties.Check(REFERENCE_HARDENING_MODULUS);
    check_properties.Check(SWELLING_SLOPE, 0.0, std::numeric_limits<double>::max());

    return 0;
}

void GeoIncrementalLinearElasticEurLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const auto& r_properties  = rValues.GetMaterialProperties();
    const auto  young_modulus = CalculateStressDependentYoungsModulus(r_properties);

    const auto poisson_ratio = r_properties.Has(POISSON_UNLOADING_RELOADING)
                                   ? r_properties[POISSON_UNLOADING_RELOADING]
                                   : r_properties[POISSON_RATIO];

    C = ConstitutiveLawUtilities::MakeContinuumElasticConstitutiveTensor(
        young_modulus, poisson_ratio, mpConstitutiveDimension->GetStrainSize(),
        mpConstitutiveDimension->GetNumberOfNormalComponents());

    if (this->GetConsiderDiagonalEntriesOnlyAndNoShear()) {
        SetEntriesAboveDiagonalToZero(C);
        SetEntriesBelowDiagonalToZero(C);
        SetShearEntriesToZero(C, mpConstitutiveDimension->GetNumberOfNormalComponents());
    }

    KRATOS_CATCH("")
}

void GeoIncrementalLinearElasticEurLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                                           Vector&       rStressVector,
                                                           ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    mDeltaStrainVector = rValues.GetStrainVector() - mStrainVectorFinalized;

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);

    noalias(mStressVector) = mStressVectorFinalized + prod(C, mDeltaStrainVector);
    rStressVector          = mStressVector;

    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticEurLaw::RequiresInitializeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticEurLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    if (!mIsModelInitialized) {
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }

    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticEurLaw::RequiresFinalizeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticEurLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void GeoIncrementalLinearElasticEurLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void GeoIncrementalLinearElasticEurLaw::ResetMaterial(const Properties&, const GeometryType&, const Vector&)
{
    mStressVector          = ZeroVector(mStressVector.size());
    mStressVectorFinalized = ZeroVector(mStressVectorFinalized.size());
    mDeltaStrainVector     = ZeroVector(mDeltaStrainVector.size());
    mStrainVectorFinalized = ZeroVector(mStrainVectorFinalized.size());
    mIsModelInitialized    = false;
}

double GeoIncrementalLinearElasticEurLaw::CalculateMinorPrincipalEffectiveStress() const
{
    auto principal_stresses = Vector{};
    auto eigen_vectors      = Matrix{};
    StressStrainUtilities::CalculatePrincipalStresses(mStressVectorFinalized, principal_stresses, eigen_vectors);

    KRATOS_DEBUG_ERROR_IF(principal_stresses.empty())
        << "Could not compute principal stresses from stress vector with size "
        << mStressVectorFinalized.size() << "\n";

    return principal_stresses[0];
}

double GeoIncrementalLinearElasticEurLaw::CalculateStressDependentYoungsModulus(const Properties& rProperties) const
{
    constexpr auto epsilon = std::numeric_limits<double>::epsilon();

    const auto reference_pressure = rProperties[REFERENCE_HARDENING_MODULUS];
    const auto exponent           = rProperties[SWELLING_SLOPE];
    const auto eur_ref            = rProperties[YOUNG_MODULUS];

    auto stress_shift = 0.0;
    if (rProperties.Has(GEO_COHESION) && rProperties.Has(GEO_FRICTION_ANGLE)) {
        const auto friction_angle_rad = MathUtils<>::DegreesToRadians(rProperties[GEO_FRICTION_ANGLE]);
        if (std::abs(std::sin(friction_angle_rad)) > epsilon) {
            stress_shift = rProperties[GEO_COHESION] * std::cos(friction_angle_rad) / std::sin(friction_angle_rad);
        }
    }

    // Keep the law stable for low confinement by enforcing -sigma3' >= p_ref.
    const auto bounded_minor_principal_effective_stress =
        std::min(CalculateMinorPrincipalEffectiveStress(), -reference_pressure);
    const auto numerator = std::max(stress_shift - bounded_minor_principal_effective_stress, epsilon);
    const auto denominator = std::max(stress_shift + reference_pressure, epsilon);

    return eur_ref * std::pow(numerator / denominator, exponent);
}

void GeoIncrementalLinearElasticEurLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
    rSerializer.save("ConstitutiveDimension", mpConstitutiveDimension);
}

void GeoIncrementalLinearElasticEurLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
    rSerializer.load("ConstitutiveDimension", mpConstitutiveDimension);
}

} // namespace Kratos
