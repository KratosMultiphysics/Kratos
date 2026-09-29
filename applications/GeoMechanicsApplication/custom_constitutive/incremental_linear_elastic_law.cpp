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
#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <string>

using namespace std::string_literals;

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

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : GeoLinearElasticLaw{},
      mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
    const auto default_properties = Properties{};
    mFormulation                  = GeoYoungsModulusFormulations::InitializeFormulation(
        GeoYoungsModulusFormulations::GetYoungsModulusFormulation(default_properties));
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
    mFormulation = rOther.mFormulation;
}

GeoIncrementalLinearElasticLaw& GeoIncrementalLinearElasticLaw::operator=(const GeoIncrementalLinearElasticLaw& rOther)
{
    if (this != &rOther) {
        GeoLinearElasticLaw::operator=(rOther);
        mStressVector          = rOther.mStressVector;
        mStressVectorFinalized = rOther.mStressVectorFinalized;
        mDeltaStrainVector     = rOther.mDeltaStrainVector;
        mStrainVectorFinalized = rOther.mStrainVectorFinalized;
        mIsModelInitialized    = rOther.mIsModelInitialized;
        if (rOther.mpConstitutiveDimension)
            mpConstitutiveDimension = rOther.mpConstitutiveDimension->Clone();
        mFormulation = rOther.mFormulation;
    }
    return *this;
}

ConstitutiveLaw::Pointer GeoIncrementalLinearElasticLaw::Clone() const
{
    return Kratos::make_shared<GeoIncrementalLinearElasticLaw>(*this);
}

bool& GeoIncrementalLinearElasticLaw::GetValue(const Variable<bool>& rVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) rValue = true;
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

int GeoIncrementalLinearElasticLaw::Check(const Properties&   rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    GeoYoungsModulusFormulations::CheckInputData(rMaterialProperties);

    return result;
}

void GeoIncrementalLinearElasticLaw::CalculateElasticMatrix(Matrix& rElasticMatrix,
                                                            ConstitutiveLaw::Parameters& rParameters)
{
    KRATOS_TRY

    const auto& r_properties = rParameters.GetMaterialProperties();

    const auto [youngs_modulus_constant, poisson_ratio] =
        ConstitutiveLawUtilities::GetOrCalculateElasticProperties(r_properties);
    const auto youngs_modulus = GeoYoungsModulusFormulations::GetYoungsModulus(
        mFormulation, r_properties, youngs_modulus_constant, mStressVectorFinalized);

    rElasticMatrix = ConstitutiveLawUtilities::MakeContinuumElasticConstitutiveTensor(
        youngs_modulus, poisson_ratio, mpConstitutiveDimension->GetStrainSize(),
        mpConstitutiveDimension->GetNumberOfNormalComponents());

    if (this->GetConsiderDiagonalEntriesOnlyAndNoShear()) {
        SetEntriesAboveDiagonalToZero(rElasticMatrix);
        SetEntriesBelowDiagonalToZero(rElasticMatrix);
        SetShearEntriesToZero(rElasticMatrix, mpConstitutiveDimension->GetNumberOfNormalComponents());
    }

    KRATOS_CATCH("")
}

void GeoIncrementalLinearElasticLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                        Vector&                      rStressVector,
                                                        ConstitutiveLaw::Parameters& rParameters)
{
    KRATOS_TRY

    mDeltaStrainVector = rParameters.GetStrainVector() - mStrainVectorFinalized;

    Matrix constitutive_matrix;
    CalculateElasticMatrix(constitutive_matrix, rParameters);

    noalias(mStressVector) = mStressVectorFinalized + prod(constitutive_matrix, mDeltaStrainVector);

    rStressVector = mStressVector;

    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticLaw::RequiresInitializeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    KRATOS_TRY
    if (!mIsModelInitialized) {
        mStressVectorFinalized = rParameters.GetStressVector();
        mStrainVectorFinalized = rParameters.GetStrainVector();
        mFormulation           = GeoYoungsModulusFormulations::InitializeFormulation(
            GeoYoungsModulusFormulations::GetYoungsModulusFormulation(rParameters.GetMaterialProperties()));
        mIsModelInitialized = true;
    }
    KRATOS_CATCH("")
}

bool GeoIncrementalLinearElasticLaw::RequiresFinalizeMaterialResponse() { return true; }

void GeoIncrementalLinearElasticLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    mStrainVectorFinalized = rParameters.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void GeoIncrementalLinearElasticLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rParameters)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rParameters);
}

void GeoIncrementalLinearElasticLaw::ResetMaterial(const Properties&, const GeometryType&, const Vector&)
{
    mStressVector          = ZeroVector(mStressVector.size());
    mStressVectorFinalized = ZeroVector(mStressVectorFinalized.size());
    mDeltaStrainVector     = ZeroVector(mDeltaStrainVector.size());
    mStrainVectorFinalized = ZeroVector(mStrainVectorFinalized.size());

    mIsModelInitialized = false;
}

void GeoIncrementalLinearElasticLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.save("ConstitutiveDimension"s, mpConstitutiveDimension);
    rSerializer.save("StressVector"s, mStressVector);
    rSerializer.save("StressVectorFinalized"s, mStressVectorFinalized);
    rSerializer.save("DeltaStrainVector"s, mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized"s, mStrainVectorFinalized);
    rSerializer.save("IsModelInitialized"s, mIsModelInitialized);
    rSerializer.save("YoungsModulusFormulations"s,
                     GeoYoungsModulusFormulations::GetYoungsModulusFormulation(mFormulation));
}

void GeoIncrementalLinearElasticLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoLinearElasticLaw)
    rSerializer.load("ConstitutiveDimension"s, mpConstitutiveDimension);
    rSerializer.load("StressVector"s, mStressVector);
    rSerializer.load("StressVectorFinalized"s, mStressVectorFinalized);
    rSerializer.load("DeltaStrainVector"s, mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized"s, mStrainVectorFinalized);
    rSerializer.load("IsModelInitialized"s, mIsModelInitialized);
    std::string youngs_modulus_formulation;
    rSerializer.load("YoungsModulusFormulations"s, youngs_modulus_formulation);
    mFormulation = GeoYoungsModulusFormulations::InitializeFormulation(youngs_modulus_formulation);
}

} // Namespace Kratos
