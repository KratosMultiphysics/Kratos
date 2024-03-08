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

// Project includes
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer GeoLinearElasticPlaneStress2DLaw::Clone() const
{
    return Kratos::make_shared<GeoLinearElasticPlaneStress2DLaw>(*this);
}

bool& GeoLinearElasticPlaneStress2DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) rValue = true;
    return rValue;
}

void GeoLinearElasticPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(PLANE_STRESS_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    // Set the strain size
    rFeatures.mStrainSize = 3;

    // Set the spacedimension
    rFeatures.mSpaceDimension = 2;
}

void GeoLinearElasticPlaneStress2DLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double      E                     = r_material_properties[YOUNG_MODULUS];
    const double      NU                    = r_material_properties[POISSON_RATIO];

    C = ZeroMatrix(GetStrainSize(), GetStrainSize());

    const double c1 = E / (1.0 - NU * NU);
    const double c2 = c1 * NU;
    const double c3 = 0.5 * E / (1.0 + NU);

    C(0, 0) = c1;
    C(0, 1) = c2;
    C(1, 0) = c2;
    C(1, 1) = c1;
    C(2, 2) = c3;
}

void GeoLinearElasticPlaneStress2DLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                                          Vector&       rStressVector,
                                                          ConstitutiveLaw::Parameters& rValues)
{
    Matrix C;
    this->CalculateElasticMatrix(C, rValues);
    noalias(rStressVector) = prod(C, rStrainVector);
}

void GeoLinearElasticPlaneStress2DLaw::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int j = 0; j < 2; ++j)
            F2x2(i, j) = F(i, j);

    Matrix E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i < 2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // namespace Kratos