// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes
// #include<cmath>

// Project includes
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicNeoHookeanPlaneStrain2D::HyperElasticIsotropicNeoHookeanPlaneStrain2D()
    : HyperElasticIsotropicNeoHookean3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticIsotropicNeoHookeanPlaneStrain2D::HyperElasticIsotropicNeoHookeanPlaneStrain2D(const HyperElasticIsotropicNeoHookeanPlaneStrain2D& rOther)
    : HyperElasticIsotropicNeoHookean3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticIsotropicNeoHookeanPlaneStrain2D::Clone() const
{
    HyperElasticIsotropicNeoHookeanPlaneStrain2D::Pointer p_clone(new HyperElasticIsotropicNeoHookeanPlaneStrain2D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicNeoHookeanPlaneStrain2D::~HyperElasticIsotropicNeoHookeanPlaneStrain2D()
{
};

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticIsotropicNeoHookeanPlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
    
    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookeanPlaneStrain2D::CalculateConstitutiveMatrixPK2(
    Matrix& ConstitutiveMatrix,
    const Matrix& InverseCTensor,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    ConstitutiveMatrix.clear();
    
    const double log_j = std::log(DeterminantF);
    
    for(unsigned int i = 0; i < 3; i++)
    {
        const unsigned int& i0 = this->msIndexVoigt2D3C[i][0];
        const unsigned int& i1 = this->msIndexVoigt2D3C[i][1];
            
        for(unsigned int j = 0; j < 3; j++)
        {
            const unsigned int& j0 = this->msIndexVoigt2D3C[j][0];
            const unsigned int& j1 = this->msIndexVoigt2D3C[j][1];
            
            ConstitutiveMatrix(i, j) = (LameLambda*InverseCTensor(i0,i1)*InverseCTensor(j0,j1)) + ((LameMu-LameLambda * log_j) * (InverseCTensor(i0,j0) * InverseCTensor(i1,j1) + InverseCTensor(i0,j1) * InverseCTensor(i1,j0)));
        }
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookeanPlaneStrain2D::CalculateConstitutiveMatrixKirchhoff(
    Matrix& ConstitutiveMatrix,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    ConstitutiveMatrix.clear();
    
    const double log_j = std::log(DeterminantF);
    
    for(unsigned int i = 0; i < 3; i++)
    {
        const unsigned int& i0 = this->msIndexVoigt2D3C[i][0];
        const unsigned int& i1 = this->msIndexVoigt2D3C[i][1];
            
        for(unsigned int j = 0; j < 3; j++)
        {
            const unsigned int& j0 = this->msIndexVoigt2D3C[j][0];
            const unsigned int& j1 = this->msIndexVoigt2D3C[j][1];
            
            ConstitutiveMatrix(i,j) = (LameLambda*((i0 == i1) ? 1.0 : 0.0)*((j0 == j1) ? 1.0 : 0.0)) + ((LameMu-LameLambda * log_j) * (((i0 == j0) ? 1.0 : 0.0) * ((i1 == j1) ? 1.0 : 0.0) + ((i0 == j1) ? 1.0 : 0.0) * ((i1 == j0) ? 1.0 : 0.0)));
        }
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookeanPlaneStrain2D::CalculateCauchyGreenStrain(
    Parameters& rValues,
    Vector& rStrainVector
    )
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // e = 0.5*(inv(C) - I)
    Matrix C_tensor = prod(trans(F),F);
    
    rStrainVector[0] = 0.5 * ( C_tensor( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( C_tensor( 1, 1 ) - 1.00 );
    rStrainVector[2] = C_tensor( 0, 1 );
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookeanPlaneStrain2D::CalculateAlmansiStrain(
    Parameters& rValues,
    Vector& rStrainVector
    )
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // e = 0.5*(1-inv(B))
    Matrix B_tensor = prod(F,trans(F));

    //Calculating the inverse of the jacobian
    Matrix inverse_B_tensor ( 2, 2 );
    double aux_det_b = 0;
    MathUtils<double>::InvertMatrix( B_tensor, inverse_B_tensor, aux_det_b);

    rStrainVector[0] = 0.5 * ( 1.0 - inverse_B_tensor( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - inverse_B_tensor( 1, 1 ) );
    rStrainVector[2] = -inverse_B_tensor( 0, 1 );
}

} // Namespace Kratos
