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
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicNeoHookean3D::HyperElasticIsotropicNeoHookean3D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticIsotropicNeoHookean3D::HyperElasticIsotropicNeoHookean3D(const HyperElasticIsotropicNeoHookean3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticIsotropicNeoHookean3D::Clone() const
{
    HyperElasticIsotropicNeoHookean3D::Pointer p_clone(new HyperElasticIsotropicNeoHookean3D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicNeoHookean3D::~HyperElasticIsotropicNeoHookean3D()
{
};

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
    
    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double& determinant_f          = rValues.GetDeterminantF();

    TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

//************************************************************************************
//************************************************************************************

void  HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();
    
    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();
    Vector& stress_vector                  = rValues.GetStressVector();
    
    // The material properties
    const double& young_modulus = material_properties[YOUNG_MODULUS];
    const double& poisson_coefficient = material_properties[POISSON_RATIO];
    
    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double& determinant_f = rValues.GetDeterminantF();
    if (determinant_f < 0.0) KRATOS_ERROR << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    
    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
    
    // Inverse of the right Cauchy-Green tensor (C):
    double aux_det;
    Matrix inverse_C_tensor(dimension, dimension); 
    MathUtils<double>::InvertMatrix( C_tensor, inverse_C_tensor, aux_det);
    
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        CalculateCauchyGreenStrain(rValues, strain_vector);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixPK2( constitutive_matrix, inverse_C_tensor, determinant_f, lame_lambda, lame_mu );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        CalculatePK2Stress( inverse_C_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();
    Vector& stress_vector                  = rValues.GetStressVector();
    
    // The material properties
    const double& young_modulus = material_properties[YOUNG_MODULUS];
    const double& poisson_coefficient = material_properties[POISSON_RATIO];
    
    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double& determinant_f = rValues.GetDeterminantF();
    if (determinant_f < 0.0) KRATOS_ERROR << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    
    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the left Cauchy-Green tensor (B):
    Matrix B_tensor = prod(deformation_gradient_f, trans( deformation_gradient_f));
    
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        CalculateAlmansiStrain(rValues, strain_vector);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixKirchoff( constitutive_matrix, determinant_f, lame_lambda, lame_mu );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        CalculateKirchoffStress( B_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);
    
    Vector& stress_vector       = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double& determinant_f = rValues.GetDeterminantF();
    
    // Set to Cauchy Stress:
    stress_vector       /= determinant_f;
    constitutive_matrix /= determinant_f;
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK1(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK2(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseCauchy(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseKirchhoff(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

double& HyperElasticIsotropicNeoHookean3D::CalculateValue(
    Parameters& rParameterValues, 
    const Variable<double>& rThisVariable, 
    double& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();
//     Vector& strain_vector                  = rParameterValues.GetStrainVector();
//     Vector& stress_vector                  = rParameterValues.GetStressVector();
    
    // The material properties
    const double& young_modulus = material_properties[YOUNG_MODULUS];
    const double& poisson_coefficient = material_properties[POISSON_RATIO];
    
    // The deformation gradient
    const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
    const double& determinant_f = rParameterValues.GetDeterminantF();
    
    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
    
    if (rThisVariable == STRAIN_ENERGY)
    {
        const double log_j = std::log(determinant_f);
        
        double first_invariant = 0.0;

        for (unsigned int i = 0; i < C_tensor.size1();i++)
        {
            first_invariant += C_tensor(i,i);
        }

        rValue = 0.5 * lame_lambda * log_j * log_j - lame_mu * log_j + 0.5 * lame_mu * (first_invariant - 3.0); 
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
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

int HyperElasticIsotropicNeoHookean3D::Check(
    const Properties& rmaterial_properties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{

    if(YOUNG_MODULUS.Key() == 0 || rmaterial_properties[YOUNG_MODULUS] <= 0.0)
    {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rmaterial_properties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
    {
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
    }

    if(DENSITY.Key() == 0 || rmaterial_properties[DENSITY] < 0.0)
    {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }

    return 0;
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateConstitutiveMatrixPK2(
    Matrix& ConstitutiveMatrix,
    const Matrix& InverseCTensor,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    ConstitutiveMatrix.clear();
    
    const double log_j = std::log(DeterminantF);
    
    for(unsigned int i = 0; i < 6; i++)
    {
        const unsigned int& i0 = this->msIndexVoigt3D6C[i][0];
        const unsigned int& i1 = this->msIndexVoigt3D6C[i][1];
            
        for(unsigned int j = 0; j < 6; j++)
        {
            const unsigned int& j0 = this->msIndexVoigt3D6C[j][0];
            const unsigned int& j1 = this->msIndexVoigt3D6C[j][1];
            
            ConstitutiveMatrix(i,j) = (LameLambda*InverseCTensor(i0,i1)*InverseCTensor(j0,j1)) + ((LameMu-LameLambda * log_j) * (InverseCTensor(i0,j0) * InverseCTensor(i1,j1) + InverseCTensor(i0,j1) * InverseCTensor(i1,j0))); 
        }
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateConstitutiveMatrixKirchoff(
    Matrix& ConstitutiveMatrix,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    ConstitutiveMatrix.clear();
    
    const double log_j = std::log(DeterminantF);
    
    for(unsigned int i = 0; i < 6; i++)
    {
        const unsigned int& i0 = this->msIndexVoigt3D6C[i][0];
        const unsigned int& i1 = this->msIndexVoigt3D6C[i][1];
            
        for(unsigned int j = 0; j < 6; j++)
        {
            const unsigned int& j0 = this->msIndexVoigt3D6C[j][0];
            const unsigned int& j1 = this->msIndexVoigt3D6C[j][1];
            
            ConstitutiveMatrix(i,j) = (LameLambda*((i0 == i1) ? 1.0 : 0.0)*((j0 == j1) ? 1.0 : 0.0)) + ((LameMu-LameLambda * log_j) * (((i0 == j0) ? 1.0 : 0.0) * ((i1 == j1) ? 1.0 : 0.0) + ((i0 == j1) ? 1.0 : 0.0) * ((i1 == j0) ? 1.0 : 0.0)));
        }
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculatePK2Stress(
    const Matrix& InvCTensor,
    Vector& rStressVector,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    Matrix stress_matrix;
    
    const SizeType dimension = WorkingSpaceDimension();
    
    stress_matrix = LameLambda * std::log(DeterminantF) * InvCTensor + LameMu * ( IdentityMatrix(dimension, dimension) - InvCTensor );
    
    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, rStressVector.size() );
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateKirchoffStress(
    const Matrix& BTensor,
    Vector& rStressVector,
    const double& DeterminantF,
    const double& LameLambda,
    const double& LameMu
    )
{
    Matrix stress_matrix;
    
    const SizeType dimension = WorkingSpaceDimension();
    
    stress_matrix  = LameLambda * std::log(DeterminantF) * IdentityMatrix(dimension, dimension) + LameMu * ( BTensor - IdentityMatrix(dimension, dimension) );
    
    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, rStressVector.size() );
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateCauchyGreenStrain(
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
    rStrainVector[2] = 0.5 * ( C_tensor( 2, 2 ) - 1.00 );
    rStrainVector[3] = C_tensor( 0, 1 ); // xy
    rStrainVector[4] = C_tensor( 1, 2 ); // yz
    rStrainVector[5] = C_tensor( 0, 2 ); // xz
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicNeoHookean3D::CalculateAlmansiStrain(
    Parameters& rValues,
    Vector& rStrainVector
    )
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // e = 0.5*(1-inv(B))
    Matrix B_tensor = prod(F,trans(F));

    //Calculating the inverse of the jacobian
    Matrix inverse_B_tensor ( 3, 3 );
    double aux_det_b = 0;
    MathUtils<double>::InvertMatrix( B_tensor, inverse_B_tensor, aux_det_b);

    rStrainVector[0] = 0.5 * ( 1.00 - inverse_B_tensor( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.00 - inverse_B_tensor( 1, 1 ) );
    rStrainVector[2] = 0.5 * ( 1.00 - inverse_B_tensor( 2, 2 ) );
    rStrainVector[3] = - inverse_B_tensor( 0, 1 ); // xy
    rStrainVector[4] = - inverse_B_tensor( 1, 2 ); // yz
    rStrainVector[5] = - inverse_B_tensor( 0, 2 ); // xz
}

} // Namespace Kratos
