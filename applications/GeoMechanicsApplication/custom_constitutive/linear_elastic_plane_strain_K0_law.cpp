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
#include "custom_constitutive/linear_elastic_plane_strain_K0_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer LinearPlaneStrainK0Law::Clone() const
{
    return Kratos::make_shared<LinearPlaneStrainK0Law>(*this);
}

bool& LinearPlaneStrainK0Law::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) rValue = true;

    return rValue;
}

void LinearPlaneStrainK0Law::GetLawFeatures(Features& rFeatures)
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

void LinearPlaneStrainK0Law::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double&     E                     = r_material_properties[YOUNG_MODULUS];
    const double&     K0ValueXX             = r_material_properties[K0_VALUE_XX];
    const double&     K0ValueYY             = r_material_properties[K0_VALUE_YY];
    const double&     K0ValueZZ             = r_material_properties[K0_VALUE_ZZ];

    double     K0_value;
    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];
    if (K0MainDirection == INDEX_2D_PLANE_STRAIN_XX) {
        K0_value = 0.5 * (K0ValueYY + K0ValueZZ);
    } else if (K0MainDirection == INDEX_2D_PLANE_STRAIN_YY) {
        K0_value = 0.5 * (K0ValueXX + K0ValueZZ);
    } else {
        KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticPlaneStrainK02DLaw: " << K0MainDirection
                     << std::endl;
    }

    double       NU    = std::max(K0_value / (K0_value + 1.0), 0.0);
    const double limit = 0.005;
    if (NU < (0.5 + limit) && NU > (0.5 - limit)) NU = 0.5 - limit;

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

void LinearPlaneStrainK0Law::CalculatePK2Stress(const Vector&                rStrainVector,
                                                Vector&                      rStressVector,
                                                ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Matrix            C;
    this->CalculateElasticMatrix(C, rValues);
    noalias(rStressVector) = prod(C, rStrainVector);

    // apply K0 procedure:
    const double& K0ValueXX = r_material_properties[K0_VALUE_XX];
    const double& K0ValueYY = r_material_properties[K0_VALUE_YY];
    const double& K0ValueZZ = r_material_properties[K0_VALUE_ZZ];

    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];
    if (K0MainDirection == INDEX_2D_PLANE_STRAIN_XX) {
        rStressVector[INDEX_2D_PLANE_STRAIN_YY] = K0ValueYY * rStressVector[INDEX_2D_PLANE_STRAIN_XX];
        rStressVector[INDEX_2D_PLANE_STRAIN_ZZ] = K0ValueZZ * rStressVector[INDEX_2D_PLANE_STRAIN_XX];
    } else if (K0MainDirection == INDEX_2D_PLANE_STRAIN_YY) {
        rStressVector[INDEX_2D_PLANE_STRAIN_XX] = K0ValueXX * rStressVector[INDEX_2D_PLANE_STRAIN_YY];
        rStressVector[INDEX_2D_PLANE_STRAIN_ZZ] = K0ValueZZ * rStressVector[INDEX_2D_PLANE_STRAIN_YY];
    } else {
        KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticPlaneStrainK02DLaw: " << K0MainDirection
                     << std::endl;
    }

    KRATOS_CATCH("")
}

void LinearPlaneStrainK0Law::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
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

} // Namespace Kratos