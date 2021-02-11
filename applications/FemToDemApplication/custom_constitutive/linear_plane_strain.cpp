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
#include "custom_constitutive/linear_plane_strain.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

LinearPlaneStrainFEMDEM::LinearPlaneStrainFEMDEM()
    : ElasticIsotropic3DFEMDEM()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

LinearPlaneStrainFEMDEM::LinearPlaneStrainFEMDEM(const LinearPlaneStrainFEMDEM& rOther)
    : ElasticIsotropic3DFEMDEM(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearPlaneStrainFEMDEM::Clone() const
{
    return Kratos::make_shared<LinearPlaneStrainFEMDEM>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

LinearPlaneStrainFEMDEM::~LinearPlaneStrainFEMDEM()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool& LinearPlaneStrainFEMDEM::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& LinearPlaneStrainFEMDEM::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& LinearPlaneStrainFEMDEM::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& LinearPlaneStrainFEMDEM::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void LinearPlaneStrainFEMDEM::GetLawFeatures(Features& rFeatures)
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

void LinearPlaneStrainFEMDEM::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
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

void LinearPlaneStrainFEMDEM::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
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

void LinearPlaneStrainFEMDEM::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i<2; ++i)
        for (unsigned int j = 0; j<2; ++j)
            F2x2(i, j) = F(i, j);

    Matrix E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i<2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos
