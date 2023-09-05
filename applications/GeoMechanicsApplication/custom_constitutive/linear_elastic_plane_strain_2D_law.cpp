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
    const double &  E = r_material_properties[YOUNG_MODULUS];
    const double & NU = r_material_properties[POISSON_RATIO];

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

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);
    // Total formulation
    noalias(rStressVector) = prod(C, rStrainVector);

    KRATOS_CATCH("")
}

} // Namespace Kratos