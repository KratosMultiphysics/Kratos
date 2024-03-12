// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes
#include <iostream>

// Project includes
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

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
    // Set the type of law
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

SizeType GeoLinearElasticPlaneStrain2DLaw::WorkingSpaceDimension() { return Dimension; }

SizeType GeoLinearElasticPlaneStrain2DLaw::GetStrainSize() const { return VoigtSize; }

bool GeoLinearElasticPlaneStrain2DLaw::IsIncremental() { return true; }

void GeoLinearElasticPlaneStrain2DLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto        E                     = r_material_properties[YOUNG_MODULUS];
    const auto        NU                    = r_material_properties[POISSON_RATIO];

    C = ZeroMatrix(GetStrainSize(), GetStrainSize());

    const double c0 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU) * c0;
    const double c2 = this->GetConsiderDiagonalEntriesOnlyAndNoShear() ? 0.0 : c0 * NU;
    const double c3 = this->GetConsiderDiagonalEntriesOnlyAndNoShear() ? 0.0 : (0.5 - NU) * c0;

    C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_XX) = c1;
    C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_YY) = c2;
    C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

    C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_XX) = c2;
    C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_YY) = c1;
    C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

    C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_XX) = c2;
    C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_YY) = c2;
    C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_ZZ) = c1;

    C(INDEX_2D_PLANE_STRAIN_XY, INDEX_2D_PLANE_STRAIN_XY) = c3;

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
        // stress and strain vectors must be initialized:
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();

        mIsModelInitialized = true;
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
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearPlaneStrainK0Law)
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("mIsModelInitialized", mIsModelInitialized);
}

void GeoLinearElasticPlaneStrain2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearPlaneStrainK0Law)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("mIsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos