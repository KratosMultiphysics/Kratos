// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Richard Faasse

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "constitutive_law_dimension.h"
#include "geo_mechanics_application_variables.h"

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

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw() = default;

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : GeoLinearElasticLaw{},
      mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw(const GeoIncrementalLinearElasticLaw& rOther)
    : GeoLinearElasticLaw(rOther),
      mStressVector(rOther.mStressVector),
      mStressVectorFinalized(rOther.mStressVectorFinalized),
      mDeltaStrainVector(rOther.mDeltaStrainVector),
      mStrainVectorFinalized(rOther.mStrainVectorFinalized),
      mIsModelInitialized(rOther.mIsModelInitialized)
{
    if (rOther.mpConstitutiveDimension)
        mpConstitutiveDimension = rOther.mpConstitutiveDimension->Clone();
}

GeoIncrementalLinearElasticLaw& GeoIncrementalLinearElasticLaw::operator=(const GeoIncrementalLinearElasticLaw& rOther)
{
    GeoLinearElasticLaw::operator=(rOther);
    mStressVector          = rOther.mStressVector;
    mStressVectorFinalized = rOther.mStressVectorFinalized;
    mDeltaStrainVector     = rOther.mDeltaStrainVector;
    mStrainVectorFinalized = rOther.mStrainVectorFinalized;
    mIsModelInitialized    = rOther.mIsModelInitialized;
    if (rOther.mpConstitutiveDimension)
        mpConstitutiveDimension = rOther.mpConstitutiveDimension->Clone();

    return *this;
}

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw(GeoIncrementalLinearElasticLaw&& rOther) noexcept = default;
GeoIncrementalLinearElasticLaw& GeoIncrementalLinearElasticLaw::operator=(GeoIncrementalLinearElasticLaw&& rOther) noexcept = default;
GeoIncrementalLinearElasticLaw::~GeoIncrementalLinearElasticLaw() = default;

ConstitutiveLaw::Pointer GeoIncrementalLinearElasticLaw::Clone() const
{
    return Kratos::make_shared<GeoIncrementalLinearElasticLaw>(*this);
}

bool& GeoIncrementalLinearElasticLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) rValue = true;
    return rValue;
}

void GeoIncrementalLinearElasticLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(mpConstitutiveDimension->GetSpatialType());
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize     = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

SizeType GeoIncrementalLinearElasticLaw::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

SizeType GeoIncrementalLinearElasticLaw::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

bool GeoIncrementalLinearElasticLaw::IsIncremental() { return true; }

void GeoIncrementalLinearElasticLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    C = mpConstitutiveDimension->CalculateElasticMatrix(r_material_properties[YOUNG_MODULUS],
                                                        r_material_properties[POISSON_RATIO]);

    if (this->GetConsiderDiagonalEntriesOnlyAndNoShear()) {
        SetEntriesAboveDiagonalToZero(C);
        SetEntriesBelowDiagonalToZero(C);
        SetShearEntriesToZero(C, mpConstitutiveDimension->GetNumberOfNormalComponents());
    }

    KRATOS_CATCH("")
}

void GeoIncrementalLinearElasticLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                        Vector&                      rStressVector,
                                                        ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    mDeltaStrainVector = rValues.GetStrainVector() - mStrainVectorFinalized;

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);

    // Incremental formulation
    noalias(mStressVector) = mStressVectorFinalized + prod(C, mDeltaStrainVector);

    rStressVector = mStressVector;

    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticLaw::RequiresInitializeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    if (!mIsModelInitialized) {
        // stress vector must be initialized:
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }
    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticLaw::RequiresFinalizeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void GeoIncrementalLinearElasticLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void GeoIncrementalLinearElasticLaw::ResetMaterial(const Properties&, const GeometryType&, const Vector&)
{
    mStressVector          = ZeroVector(mStressVector.size());
    mStressVectorFinalized = ZeroVector(mStressVectorFinalized.size());

    // set strain vectors:
    mDeltaStrainVector     = ZeroVector(mDeltaStrainVector.size());
    mStrainVectorFinalized = ZeroVector(mStrainVectorFinalized.size());

    mIsModelInitialized = false;
}

void GeoIncrementalLinearElasticLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("mIsModelInitialized", mIsModelInitialized);
}

void GeoIncrementalLinearElasticLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("mIsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos
