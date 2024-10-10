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

#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"
#include "constitutive_law_dimension.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

GeoLinearElasticPlaneStrain2DLaw::GeoLinearElasticPlaneStrain2DLaw() = default;

GeoLinearElasticPlaneStrain2DLaw::GeoLinearElasticPlaneStrain2DLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : GeoLinearElasticLaw{},
      mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

GeoLinearElasticPlaneStrain2DLaw::GeoLinearElasticPlaneStrain2DLaw(const GeoLinearElasticPlaneStrain2DLaw& rOther)
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

GeoLinearElasticPlaneStrain2DLaw& GeoLinearElasticPlaneStrain2DLaw::operator=(const GeoLinearElasticPlaneStrain2DLaw& rOther)
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

GeoLinearElasticPlaneStrain2DLaw::GeoLinearElasticPlaneStrain2DLaw(GeoLinearElasticPlaneStrain2DLaw&& rOther) = default;
GeoLinearElasticPlaneStrain2DLaw& GeoLinearElasticPlaneStrain2DLaw::operator=(GeoLinearElasticPlaneStrain2DLaw&& rOther) = default;
GeoLinearElasticPlaneStrain2DLaw::~GeoLinearElasticPlaneStrain2DLaw() = default;

ConstitutiveLaw::Pointer GeoLinearElasticPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<GeoLinearElasticPlaneStrain2DLaw>(*this);
}

bool& GeoLinearElasticPlaneStrain2DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) rValue = true;
    return rValue;
}

void GeoLinearElasticPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(mpConstitutiveDimension->GetSpatialType());
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize     = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

SizeType GeoLinearElasticPlaneStrain2DLaw::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

SizeType GeoLinearElasticPlaneStrain2DLaw::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

bool GeoLinearElasticPlaneStrain2DLaw::IsIncremental() { return true; }

void GeoLinearElasticPlaneStrain2DLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto        E                     = r_material_properties[YOUNG_MODULUS];
    const auto        NU                    = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU) * c0;
    const double c2 = this->GetConsiderDiagonalEntriesOnlyAndNoShear() ? 0.0 : c0 * NU;
    const double c3 = this->GetConsiderDiagonalEntriesOnlyAndNoShear() ? 0.0 : (0.5 - NU) * c0;

    C = mpConstitutiveDimension->FillConstitutiveMatrix(c1, c2, c3);

    KRATOS_CATCH("")
}

void GeoLinearElasticPlaneStrain2DLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                                          Vector&       rStressVector,
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

bool GeoLinearElasticPlaneStrain2DLaw::RequiresInitializeMaterialResponse() { return true; }

void GeoLinearElasticPlaneStrain2DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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

bool GeoLinearElasticPlaneStrain2DLaw::RequiresFinalizeMaterialResponse() { return true; }

void GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void GeoLinearElasticPlaneStrain2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("mIsModelInitialized", mIsModelInitialized);
}

void GeoLinearElasticPlaneStrain2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("mIsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos