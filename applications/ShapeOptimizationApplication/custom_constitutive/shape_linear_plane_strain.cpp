// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/shape_linear_plane_strain.h"

#include "shape_optimization_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

ShapeLinearPlaneStrain::ShapeLinearPlaneStrain()
    : ShapeElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

ShapeLinearPlaneStrain::ShapeLinearPlaneStrain(const ShapeLinearPlaneStrain& rOther)
    : ShapeElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer ShapeLinearPlaneStrain::Clone() const
{
    return Kratos::make_shared<ShapeLinearPlaneStrain>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

ShapeLinearPlaneStrain::~ShapeLinearPlaneStrain()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool& ShapeLinearPlaneStrain::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == SHAPE_STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ShapeLinearPlaneStrain::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

// ConstitutiveLaw::VoigtSizeMatrixType& ShapeLinearPlaneStrain::GetValue(const Variable<VoigtSizeMatrixType>& rThisVariable, ConstitutiveLaw::VoigtSizeMatrixType& rValue)
// {
//     return rValue;
// }

/***********************************************************************************/
/***********************************************************************************/

// ConstitutiveLaw::DeformationGradientMatrixType& ShapeLinearPlaneStrain::GetValue(const Variable<DeformationGradientMatrixType>& rThisVariable, ConstitutiveLaw::DeformationGradientMatrixType& rValue)
// {
//     return rValue;
// }

/***********************************************************************************/
/***********************************************************************************/

Vector& ShapeLinearPlaneStrain::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

// ConstitutiveLaw::StrainVectorType& ShapeLinearPlaneStrain::GetValue(const Variable<StrainVectorType>& rThisVariable, ConstitutiveLaw::StrainVectorType& rValue)
// {
//     return rValue;
// }

/***********************************************************************************/
/***********************************************************************************/

// ConstitutiveLaw::StressVectorType& ShapeLinearPlaneStrain::GetValue(const Variable<StressVectorType>& rThisVariable, ConstitutiveLaw::StressVectorType& rValue)
// {
//     return rValue;
// }

/***********************************************************************************/
/***********************************************************************************/

double& ShapeLinearPlaneStrain::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void ShapeLinearPlaneStrain::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 3;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 2;
}

/***********************************************************************************/
/***********************************************************************************/

void ShapeLinearPlaneStrain::CalculateElasticMatrix(VoigtSizeMatrixType& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(C);

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
    const double c1 = (1.00 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    C(0, 0) = c1;
    C(0, 1) = c2;
    C(1, 0) = c2;
    C(1, 1) = c1;
    C(2, 2) = c3;
}

/***********************************************************************************/
/***********************************************************************************/

void ShapeLinearPlaneStrain::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
    const double c1 = (1.00 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    rStressVector[0] = c1 * rStrainVector[0] + c2 * rStrainVector[1];
    rStressVector[1] = c2 * rStrainVector[0] + c1 * rStrainVector[1];
    rStressVector[2] = c3 * rStrainVector[2];
}

/***********************************************************************************/
/***********************************************************************************/

void ShapeLinearPlaneStrain::CalculateCauchyGreenStrain(Parameters& rValues, ConstitutiveLaw::StrainVectorType& rStrainVector)
{
    //1.-Compute total deformation gradient
    const ConstitutiveLaw::DeformationGradientMatrixType& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i<2; ++i)
        for (unsigned int j = 0; j<2; ++j)
            F2x2(i, j) = F(i, j);

    BoundedMatrix<double,2,2> E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i<2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos
