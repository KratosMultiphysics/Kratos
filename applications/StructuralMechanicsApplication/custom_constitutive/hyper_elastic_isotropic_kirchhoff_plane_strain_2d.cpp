// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Malik Ali Dawi
//                   Ruben Zorrilla
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
HyperElasticIsotropicKirchhoffPlaneStrain2D::HyperElasticIsotropicKirchhoffPlaneStrain2D()
    : ConstitutiveLaw() {
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticIsotropicKirchhoffPlaneStrain2D::HyperElasticIsotropicKirchhoffPlaneStrain2D(const HyperElasticIsotropicKirchhoffPlaneStrain2D& rOther)
    : ConstitutiveLaw(rOther) {
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticIsotropicKirchhoffPlaneStrain2D::Clone() const {
    HyperElasticIsotropicKirchhoffPlaneStrain2D::Pointer p_clone(new HyperElasticIsotropicKirchhoffPlaneStrain2D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicKirchhoffPlaneStrain2D::~HyperElasticIsotropicKirchhoffPlaneStrain2D() {
};

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) {

    CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double& determinant_f          = rValues.GetDeterminantF();

    TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

//************************************************************************************
//************************************************************************************

void  HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    // Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();

    // The material properties
    const double& young_modulus = material_properties[YOUNG_MODULUS];
    const double& poisson_coefficient = material_properties[POISSON_RATIO];


    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateGreenLagrangianStrain(rValues, strain_vector);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixPK2( ConstitutiveMatrix, young_modulus, poisson_coefficient);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector& stress_vector = rValues.GetStressVector();
        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
            noalias(stress_vector) = prod(constitutive_matrix, strain_vector);
        } else {
            CalculatePK2Stress( strain_vector, stress_vector, young_modulus, poisson_coefficient );
        }
    }
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) {
    
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

    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateAlmansiStrain(rValues, strain_vector);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixKirchhoff(constitutive_matrix, deformation_gradient_f ,young_modulus, poisson_coefficient);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        if (rValues.IsSetDeformationGradientF() == true) {
               CalculateGreenLagrangianStrain(rValues, strain_vector);
        }
        CalculateKirchhoffStress(strain_vector, stress_vector, deformation_gradient_f ,young_modulus, poisson_coefficient);
    }
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) {

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

void HyperElasticIsotropicKirchhoffPlaneStrain2D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK1(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK2(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseCauchy(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseKirchhoff(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

//************************************************************************************
//************************************************************************************


double& HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue) {

    const Properties& material_properties  = rParameterValues.GetMaterialProperties();
    Vector& strain_vector                  = rParameterValues.GetStrainVector();

    // The material properties
    const double& young_modulus = material_properties[YOUNG_MODULUS];
    const double& poisson_coefficient = material_properties[POISSON_RATIO];

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    if (rThisVariable == STRAIN_ENERGY) {
        CalculateGreenLagrangianStrain(rParameterValues,strain_vector);
        Matrix E_tensor=MathUtils<double>::StrainVectorToTensor(strain_vector);
        Matrix E_tensor_sq=prod(E_tensor,E_tensor);
        double E_trace = 0.0;
        double E_trace_sq = 0.0;
        for (unsigned int i = 0; i < E_tensor.size1();i++) {
            E_trace     += E_tensor (i,i);
            E_trace_sq  += E_tensor_sq(i,i);
        }

        rValue = 0.5 * lame_lambda * E_trace * E_trace  + 0.5 *lame_mu *E_trace_sq ;
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::GetLawFeatures(Features& rFeatures) {

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

int HyperElasticIsotropicKirchhoffPlaneStrain2D::Check(
    const Properties& rmaterial_properties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) {

    if(YOUNG_MODULUS.Key() == 0 || rmaterial_properties[YOUNG_MODULUS] <= 0.0) {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rmaterial_properties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true) {
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
    }

    if(DENSITY.Key() == 0 || rmaterial_properties[DENSITY] < 0.0) {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }

    return 0;
}

//************************************************************************************
//************************************************************************************
//CalculateConstitutiveMatrixPK2( constitutive_matrix, inverse_C_tensor, determinant_f, lame_lambda, lame_mu );

//CalculateConstitutiveMatrixPK2( ConstitutiveMatrix, young_modulus, poisson_coefficient);
void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateConstitutiveMatrixPK2(
    Matrix& ConstitutiveMatrix,
    const double& YoungModulus,
    const double& PoissonCoefficient) {

    ConstitutiveMatrix.clear();
    ConstitutiveMatrix = ZeroMatrix(3,3);
    const double c0 = YoungModulus / ((1.00 + PoissonCoefficient)*(1-2*PoissonCoefficient));
    const double c1 = (1.00 - PoissonCoefficient)*c0;
    const double c2 = c0 * PoissonCoefficient;
    const double c3 = (0.5 - PoissonCoefficient)*c0;

    ConstitutiveMatrix(0,0) = c1;
    ConstitutiveMatrix(0,1) = c2;
    ConstitutiveMatrix(1,0) = c2;
    ConstitutiveMatrix(1,1) = c1;
    ConstitutiveMatrix(2,2) = c3;
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateConstitutiveMatrixKirchhoff(
    Matrix& ConstitutiveMatrix,
    const Matrix& DeformationGradientF,
    const double& YoungModulus,
    const double& PoissonCoefficient) {

    ConstitutiveMatrix.clear();

    CalculateConstitutiveMatrixPK2(ConstitutiveMatrix, YoungModulus, PoissonCoefficient);
    PushForwardConstitutiveMatrix (ConstitutiveMatrix, DeformationGradientF );
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    const double& YoungModulus,
    const double& PoissonCoefficient) {

    SizeType SizeSystem = GetStrainSize();
    Matrix ConstMatrixPK2(SizeSystem,SizeSystem);
    CalculateConstitutiveMatrixPK2(ConstMatrixPK2, YoungModulus, PoissonCoefficient);
    rStressVector = prod(ConstMatrixPK2,rStrainVector);
}


//************************************************************************************
//************************************************************************************
void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateKirchhoffStress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    const Matrix& DeformationGradientF,
    const double& YoungModulus,
    const double& PoissonCoefficient) {

    CalculatePK2Stress( rStrainVector, rStressVector, YoungModulus, PoissonCoefficient );
    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
    ContraVariantPushForward (StressMatrix,DeformationGradientF); //Kirchhoff
    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector) {

    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // E = 0.5*(inv(C) - I)
    Matrix C_tensor = prod(trans(F),F);

    rStrainVector[0] = 0.5 * ( C_tensor( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( C_tensor( 1, 1 ) - 1.00 );
    rStrainVector[2] = C_tensor( 0, 1 ); // xy
}

//************************************************************************************
//************************************************************************************

void HyperElasticIsotropicKirchhoffPlaneStrain2D::CalculateAlmansiStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector) {
        
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // e = 0.5*(1-inv(B))
    Matrix B_tensor = prod(F,trans(F));

    //Calculating the inverse of the jacobian
    Matrix inverse_B_tensor ( 2, 2 );
    double aux_det_b = 0;
    MathUtils<double>::InvertMatrix( B_tensor, inverse_B_tensor, aux_det_b);

    rStrainVector[0] = 0.5 * ( 1.00 - inverse_B_tensor( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.00 - inverse_B_tensor( 1, 1 ) );
    rStrainVector[2] = - inverse_B_tensor( 0, 1 ); // xy ??? yz
}

} // Namespace Kratos
