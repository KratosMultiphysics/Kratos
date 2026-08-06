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
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>
#include <limits>
#include <string>

using namespace std::string_literals;

namespace Kratos
{

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

GeoIncrementalLinearElasticLaw::GeoIncrementalLinearElasticLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : GeoLinearElasticLaw{},
      mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
    mPolicy.emplace<Policies::Constant>();
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
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    if (rMaterialProperties.Has(GEO_YOUNGS_MODULUS_FORMULATION)) {
        if (rMaterialProperties[GEO_YOUNGS_MODULUS_FORMULATION] == "Eur") {
            const CheckProperties check_properties(rMaterialProperties, "parameters of material",
                                                   CheckProperties::Bounds::AllExclusive);
            check_properties.Check(GEO_PRESSURE_REFERENCE);
            check_properties.Check(GEO_STRESS_DEPENDENCY_EXPONENT);
            check_properties.Check(GEO_COHESION);
            check_properties.Check(GEO_FRICTION_ANGLE);
        }
    }

    return 0;
}

void GeoIncrementalLinearElasticLaw::CalculateElasticMatrix(Matrix& rElasticMatrix,
                                                            ConstitutiveLaw::Parameters& rParameters)
{
    KRATOS_TRY

    const auto& r_properties   = rParameters.GetMaterialProperties();
    const auto  youngs_modulus = GetYoungsModulus(r_properties);

    const auto poisson_ratio = r_properties.Has(POISSON_UNLOADING_RELOADING)
                                   ? r_properties[POISSON_UNLOADING_RELOADING]
                                   : r_properties[POISSON_RATIO];

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

double GeoIncrementalLinearElasticLaw::GetYoungsModulus(const Properties& rProperties) const
{
    return std::visit([this, &rProperties](auto&& policy) -> double {
        using T = std::decay_t<decltype(policy)>;

        if constexpr (std::is_same_v<T, Policies::Constant>) {
            return rProperties[YOUNG_MODULUS];

        } else if constexpr (std::is_same_v<T, Policies::Eur>) {
            return CalculateYoungsModulusForEur(rProperties);
        }
    }, mPolicy);
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

        const auto& r_material_props = rParameters.GetMaterialProperties();
        InitializePolicy(r_material_props);

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

    // set strain vectors:
    mDeltaStrainVector     = ZeroVector(mDeltaStrainVector.size());
    mStrainVectorFinalized = ZeroVector(mStrainVectorFinalized.size());

    mIsModelInitialized = false;
}

void GeoIncrementalLinearElasticLaw::InitializePolicy(const Properties& rProperties)
{
    if (rProperties.Has(GEO_YOUNGS_MODULUS_FORMULATION)) {
        const std::string& formulation = rProperties[GEO_YOUNGS_MODULUS_FORMULATION];
        if (formulation == "Constant") {
            mPolicy.emplace<Policies::Constant>();
        } else if (formulation == "Eur") {
            mPolicy.emplace<Policies::Eur>();
        } else {
            KRATOS_ERROR << "Unknown GEO_YOUNGS_MODULUS_FORMULATION: " << formulation;
        }
    } else {
        // Default fallback
        mPolicy.emplace<Policies::Constant>();
    }
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
}

double GeoIncrementalLinearElasticLaw::CalculateYoungsModulusForEur(const Properties& rProperties) const
{
    constexpr auto epsilon = std::numeric_limits<double>::epsilon();

    const auto reference_pressure = rProperties[GEO_PRESSURE_REFERENCE];
    const auto exponent           = rProperties[GEO_STRESS_DEPENDENCY_EXPONENT];
    const auto eur_ref            = rProperties[YOUNG_MODULUS];

    const auto friction_angle_rad = ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties);
    const auto stress_shift =
        rProperties[GEO_COHESION] * std::cos(friction_angle_rad) / std::sin(friction_angle_rad);

    const auto base =
        (stress_shift - CalculateMinorPrincipalEffectiveStress()) / (stress_shift + reference_pressure);

    KRATOS_ERROR_IF_NOT(base > epsilon)
        << "Negative base for std::pow ("
        << base << "). Check GEO_COHESION, GEO_FRICTION_ANGLE, GEO_PRESSURE_REFERENCE and the finalized stress state.\n";

    return eur_ref * std::pow(base, exponent);
}

double GeoIncrementalLinearElasticLaw::CalculateMinorPrincipalEffectiveStress() const
{
    auto principal_stresses = Vector{};
    auto eigen_vectors      = Matrix{};
    StressStrainUtilities::CalculatePrincipalStresses(mStressVectorFinalized, principal_stresses, eigen_vectors);

    KRATOS_DEBUG_ERROR_IF(principal_stresses.size() < 3)
        << "Could not compute principal stresses from stress vector with size "
        << mStressVectorFinalized.size() << ". Expected at least 3 principal stresses, got "
        << principal_stresses.size() << "\n";

    return principal_stresses[2];
}

} // Namespace Kratos
