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
#include "custom_constitutive/linear_plane_stress_pre_strain.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStressPreStrain::LinearPlaneStressPreStrain()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearPlaneStressPreStrain::LinearPlaneStressPreStrain(const LinearPlaneStressPreStrain& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearPlaneStressPreStrain::Clone() const
{
    LinearPlaneStressPreStrain::Pointer p_clone(new LinearPlaneStressPreStrain(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStressPreStrain::~LinearPlaneStressPreStrain()
{
}

//************************************************************************************
//************************************************************************************

bool& LinearPlaneStressPreStrain::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE)
        rValue = true;

    return rValue;
}


bool LinearPlaneStressPreStrain::Has(const Variable<double>& rThisVariable)

{
    if (rThisVariable == PRE_STRAIN_FACTOR) {
        return true;
    }
    return false;
}

void LinearPlaneStressPreStrain::SetValue(const Variable<double>& rVariable,
    const double& Value,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == PRE_STRAIN_FACTOR) {
        m_pre_strain_factor = Value;
    }
}
//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearPlaneStressPreStrain::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
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

//************************************************************************************
//************************************************************************************

void LinearPlaneStressPreStrain::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(C);

    const double c1 = E / (1.00 - NU * NU);
    const double c2 = c1 * NU;
    const double c3 = 0.5*E / (1 + NU);

    C(0, 0) = c1;
    C(0, 1) = c2;
    C(1, 0) = c2;
    C(1, 1) = c1;
    C(2, 2) = c3;
}

//************************************************************************************
//************************************************************************************

void LinearPlaneStressPreStrain::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c1 = E / (1.00 - NU * NU);
    const double c2 = c1 * NU;
    const double c3 = 0.5* E / (1 + NU);

    Vector strain_with_pre_strain = rStrainVector;
    strain_with_pre_strain[0] += m_pre_strain_factor;
    strain_with_pre_strain[1] += m_pre_strain_factor;

    rStressVector[0] = c1 * strain_with_pre_strain[0] + c2 * strain_with_pre_strain[1];
    rStressVector[1] = c2 * strain_with_pre_strain[0] + c1 * strain_with_pre_strain[1];
    rStressVector[2] = c3 * strain_with_pre_strain[2];
}

//************************************************************************************
//************************************************************************************

void LinearPlaneStressPreStrain::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
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
