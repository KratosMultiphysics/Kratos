// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//

// Project includes
#include "custom_constitutive/elastic_isotropic_K0_3d_law.h"
#include "includes/checks.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer ElasticIsotropicK03DLaw::Clone() const
{
    return Kratos::make_shared<ElasticIsotropicK03DLaw>(*this);
}

void ElasticIsotropicK03DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
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

void ElasticIsotropicK03DLaw::CalculateElasticMatrix(Matrix& rConstitutiveMatrix, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double      E                     = r_material_properties[YOUNG_MODULUS];

    const double& K0ValueXX = r_material_properties[K0_VALUE_XX];
    const double& K0ValueYY = r_material_properties[K0_VALUE_YY];
    const double& K0ValueZZ = r_material_properties[K0_VALUE_ZZ];

    double     K0_value;
    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];
    if (K0MainDirection == INDEX_3D_XX) {
        K0_value = 0.5 * (K0ValueYY + K0ValueZZ);
    } else if (K0MainDirection == INDEX_3D_YY) {
        K0_value = 0.5 * (K0ValueXX + K0ValueZZ);
    } else if (K0MainDirection == INDEX_3D_ZZ) {
        K0_value = 0.5 * (K0ValueXX + K0ValueYY);
    } else {
        KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticK03DLaw: " << K0MainDirection << std::endl;
    }

    double       NU    = std::max(K0_value / (K0_value + 1.0), 0.0);
    const double limit = 0.005;
    if (NU < (0.5 + limit) && NU > (0.5 - limit)) NU = 0.5 - limit;

    const double c1 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c2 = c1 * (1.0 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1.0 - 2.0 * NU);

    rConstitutiveMatrix       = ZeroMatrix(GetStrainSize(), GetStrainSize());
    rConstitutiveMatrix(0, 0) = c2;
    rConstitutiveMatrix(0, 1) = c3;
    rConstitutiveMatrix(0, 2) = c3;
    rConstitutiveMatrix(1, 0) = c3;
    rConstitutiveMatrix(1, 1) = c2;
    rConstitutiveMatrix(1, 2) = c3;
    rConstitutiveMatrix(2, 0) = c3;
    rConstitutiveMatrix(2, 1) = c3;
    rConstitutiveMatrix(2, 2) = c2;
    rConstitutiveMatrix(3, 3) = c4;
    rConstitutiveMatrix(4, 4) = c4;
    rConstitutiveMatrix(5, 5) = c4;
}

void ElasticIsotropicK03DLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                 Vector&                      rStressVector,
                                                 ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Matrix            C;
    this->CalculateElasticMatrix(C, rValues);
    noalias(rStressVector) = prod(C, rStrainVector);

    // apply K0 procedure:
    const double& K0ValueXX = r_material_properties[K0_VALUE_XX];
    const double& K0ValueYY = r_material_properties[K0_VALUE_YY];
    const double& K0ValueZZ = r_material_properties[K0_VALUE_ZZ];

    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];
    if (K0MainDirection == INDEX_3D_XX) {
        rStressVector[INDEX_3D_YY] = K0ValueYY * rStressVector[INDEX_3D_XX];
        rStressVector[INDEX_3D_ZZ] = K0ValueZZ * rStressVector[INDEX_3D_XX];
    } else if (K0MainDirection == INDEX_3D_YY) {
        rStressVector[INDEX_3D_XX] = K0ValueXX * rStressVector[INDEX_3D_YY];
        rStressVector[INDEX_3D_ZZ] = K0ValueZZ * rStressVector[INDEX_3D_YY];
    } else if (K0MainDirection == INDEX_3D_ZZ) {
        rStressVector[INDEX_3D_XX] = K0ValueXX * rStressVector[INDEX_3D_ZZ];
        rStressVector[INDEX_3D_YY] = K0ValueYY * rStressVector[INDEX_3D_ZZ];
    } else {
        KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticK03DLaw: " << K0MainDirection << std::endl;
    }
}

void ElasticIsotropicK03DLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector)
{
    const SizeType space_dimension = this->WorkingSpaceDimension();

    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1() != space_dimension || F.size2() != space_dimension)
        << "expected size of F " << space_dimension << "x" << space_dimension << ", got "
        << F.size1() << "x" << F.size2() << std::endl;

    Matrix E_tensor = prod(trans(F), F);
    for (unsigned int i = 0; i < space_dimension; ++i)
        E_tensor(i, i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos