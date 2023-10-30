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
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

void GeoLinearElasticPlaneStrain2DLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties &r_material_properties = rValues.GetMaterialProperties();
    const double &  E  = r_material_properties[YOUNG_MODULUS];
    const double &  NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(C);

    const double c0 = E / ((1.0 + NU)*(1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

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
                                                          Vector& rStressVector,
                                                          ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    UpdateInternalDeltaStrainVector(rValues);

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);

    // Incremental formulation
    noalias(mStressVector) = mStressVectorFinalized + prod(C, mDeltaStrainVector);
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::CalculatePK2Stress") <<
         "mDeltaStrainVector = " << mDeltaStrainVector[0] << ", " << mDeltaStrainVector[1] << ", " <<  mDeltaStrainVector[2] << ", " <<  mDeltaStrainVector[3] << std::endl;
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::CalculatePK2Stress") <<
         "mStressVectorFinalized = " << mStressVectorFinalized[0] << ", " << mStressVectorFinalized[1] << ", " <<  mStressVectorFinalized[2] << ", " <<  mStressVectorFinalized[3] << std::endl;
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::CalculatePK2Stress") <<
         "mStressVector = " << mStressVector[0] << ", " << mStressVector[1] << ", " <<  mStressVector[2] << ", " <<  mStressVector[3] << std::endl;

    SetExternalStressVector(rStressVector);

    KRATOS_CATCH("")
}

void GeoLinearElasticPlaneStrain2DLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues){
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::UpdateInternalDeltaStrainVector") << std::endl;
    const Vector& rStrainVector = rValues.GetStrainVector();
    for (unsigned int i=0; i < mDeltaStrainVector.size(); ++i)
    {
        mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
    }
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::UpdateInternalDeltaStrainVector") <<
        "mDeltaStrainVector = " << mDeltaStrainVector[0] << ", " << mDeltaStrainVector[1] << ", " <<  mDeltaStrainVector[2] << ", " <<  mDeltaStrainVector[3] << std::endl;
}

void GeoLinearElasticPlaneStrain2DLaw::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters &rValues) {
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::UpdateInternalStrainVectorFinalized") << std::endl;
    const Vector& rStrainVector = rValues.GetStrainVector();
    this->SetInternalStrainVector(rStrainVector);
}

void GeoLinearElasticPlaneStrain2DLaw::SetExternalStressVector(Vector& rStressVector){
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::SetExternalStressVector") << std::endl;
    for (unsigned int i=0; i < rStressVector.size(); ++i)
    {
        rStressVector(i) = mStressVector[i];
    }
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::SetExternalStressVector") <<
      "mStressVector = " << mStressVector[0] << ", " << mStressVector[1] << ", " <<  mStressVector[2] << ", " <<  mStressVector[3] << std::endl;
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::SetExternalStressVector") <<
       "rStressVector = " << rStressVector[0] << ", " << rStressVector[1] << ", " <<  rStressVector[2] << ", " <<  rStressVector[3] << std::endl;
}

void GeoLinearElasticPlaneStrain2DLaw::SetInternalStressVector(const Vector& rStressVector){
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::SetInternalStressVector") << std::endl;
    for (unsigned int i=0; i < mStressVectorFinalized.size(); ++i)
    {
        mStressVectorFinalized[i] = rStressVector(i);
    }
}

void GeoLinearElasticPlaneStrain2DLaw::SetInternalStrainVector(const Vector& rStrainVector){
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::SetInternalStrainVector") << std::endl;
    for (unsigned int i=0; i < mStrainVectorFinalized.size(); ++i)
    {
        mStrainVectorFinalized[i] = rStrainVector(i);
    }
}
void GeoLinearElasticPlaneStrain2DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::InitializeMaterialResponseCauchy") << std::endl;
    KRATOS_TRY
    if (!mIsModelInitialized)
    {
        // stress and strain vectors must be initialized:
        const Vector& rStressVector = rValues.GetStressVector();
        const Vector& rStrainVector = rValues.GetStrainVector();

        SetInternalStressVector(rStressVector);
        SetInternalStrainVector(rStrainVector);

        mIsModelInitialized = true;
    }
    KRATOS_CATCH("")
}

void GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters & rValues)
{
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponseCauchy") << std::endl;
    UpdateInternalStrainVectorFinalized(rValues);
    mStressVectorFinalized = mStressVector;
}

void GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    KRATOS_INFO("GeoLinearElasticPlaneStrain2DLaw::FinalizeMaterialResponsePK2") << std::endl;
    FinalizeMaterialResponseCauchy(rValues);
}


} // Namespace Kratos