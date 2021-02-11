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

// Project includes
#include "includes/checks.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicNeoHookean3D::HyperElasticIsotropicNeoHookean3D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicNeoHookean3D::HyperElasticIsotropicNeoHookean3D(const HyperElasticIsotropicNeoHookean3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicNeoHookean3D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicNeoHookean3D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicNeoHookean3D::~HyperElasticIsotropicNeoHookean3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void  HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);

    // Inverse of the right Cauchy-Green tensor (C):
    double aux_det;
    Matrix inverse_C_tensor(dimension, dimension);
    MathUtils<double>::InvertMatrix( C_tensor, inverse_C_tensor, aux_det);

    if(r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateGreenLagrangianStrain(rValues, strain_vector);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixPK2( constitutive_matrix, inverse_C_tensor, determinant_f, lame_lambda, lame_mu );
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector& stress_vector = rValues.GetStressVector();
        CalculatePK2Stress( inverse_C_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();
    Vector& stress_vector                  = rValues.GetStressVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    if(r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateAlmansiStrain(rValues, strain_vector);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixKirchhoff( constitutive_matrix, determinant_f, lame_lambda, lame_mu );
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        // We compute the left Cauchy-Green tensor (B):
        const Matrix B_tensor = prod(deformation_gradient_f, trans( deformation_gradient_f));
        CalculateKirchhoffStress( B_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);

    Vector& stress_vector       = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f  = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    stress_vector       /= determinant_f;
    constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

double& HyperElasticIsotropicNeoHookean3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
    const double determinant_f = rParameterValues.GetDeterminantF();

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);

    if (rThisVariable == STRAIN_ENERGY) {
        const double log_j = std::log(determinant_f);

        double first_invariant = 0.0;

        for (IndexType i = 0; i < C_tensor.size1();++i) {
            first_invariant += C_tensor(i,i);
        }

        rValue = 0.5 * lame_lambda * log_j * log_j - lame_mu * log_j + 0.5 * lame_mu * (first_invariant - 3.0);
    }

    return( rValue );
}


/***********************************************************************************/
/***********************************************************************************/

Vector& HyperElasticIsotropicNeoHookean3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == HENCKY_STRAIN_VECTOR ||
        rThisVariable == BIOT_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the strain
        if (rThisVariable == STRAIN) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == ALMANSI_STRAIN_VECTOR) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseKirchhoff(rParameterValues);
        } else if (rThisVariable == HENCKY_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_tensor, r_strain_vector);
        } else if (rThisVariable == BIOT_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateBiotStrain(C_tensor, r_strain_vector);
        }

        rValue = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        if (rThisVariable == STRESSES) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseKirchhoff(rParameterValues);
        } if (rThisVariable == CAUCHY_STRESS_VECTOR) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponseCauchy(rParameterValues);
        } if (rThisVariable == PK2_STRESS_VECTOR) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& HyperElasticIsotropicNeoHookean3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            HyperElasticIsotropicNeoHookean3D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

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

/***********************************************************************************/
/***********************************************************************************/

int HyperElasticIsotropicNeoHookean3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is null or negative." << std::endl;

    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateConstitutiveMatrixPK2(
    Matrix& rConstitutiveMatrix,
    const Matrix& rInverseCTensor,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    rConstitutiveMatrix.clear();

    const double log_j = std::log(DeterminantF);

    for(IndexType i = 0; i < 6; ++i) {
        const IndexType i0 = this->msIndexVoigt3D6C[i][0];
        const IndexType i1 = this->msIndexVoigt3D6C[i][1];

        for(IndexType j = 0; j < 6; ++j) {
            const IndexType j0 = this->msIndexVoigt3D6C[j][0];
            const IndexType j1 = this->msIndexVoigt3D6C[j][1];

            rConstitutiveMatrix(i,j) = (LameLambda*rInverseCTensor(i0,i1)*rInverseCTensor(j0,j1)) + ((LameMu-LameLambda * log_j) * (rInverseCTensor(i0,j0) * rInverseCTensor(i1,j1) + rInverseCTensor(i0,j1) * rInverseCTensor(i1,j0)));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateConstitutiveMatrixKirchhoff(
    Matrix& rConstitutiveMatrix,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    rConstitutiveMatrix.clear();

    const double log_j = std::log(DeterminantF);

    for(IndexType i = 0; i < 6; ++i) {
        const IndexType i0 = this->msIndexVoigt3D6C[i][0];
        const IndexType i1 = this->msIndexVoigt3D6C[i][1];

        for(IndexType j = 0; j < 6; ++j) {
            const IndexType j0 = this->msIndexVoigt3D6C[j][0];
            const IndexType j1 = this->msIndexVoigt3D6C[j][1];

            rConstitutiveMatrix(i,j) = (LameLambda*((i0 == i1) ? 1.0 : 0.0)*((j0 == j1) ? 1.0 : 0.0)) + ((LameMu-LameLambda * log_j) * (((i0 == j0) ? 1.0 : 0.0) * ((i1 == j1) ? 1.0 : 0.0) + ((i0 == j1) ? 1.0 : 0.0) * ((i1 == j0) ? 1.0 : 0.0)));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculatePK2Stress(
    const Matrix& rInvCTensor,
    Vector& rStressVector,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    Matrix stress_matrix;

    const SizeType dimension = WorkingSpaceDimension();

    stress_matrix = LameLambda * std::log(DeterminantF) * rInvCTensor + LameMu * ( IdentityMatrix(dimension) - rInvCTensor );

    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, GetStrainSize() );
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateKirchhoffStress(
    const Matrix& rBTensor,
    Vector& rStressVector,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    Matrix stress_matrix;

    const SizeType dimension = WorkingSpaceDimension();

    stress_matrix  = LameLambda * std::log(DeterminantF) * IdentityMatrix(dimension) + LameMu * ( rBTensor - IdentityMatrix(dimension) );

    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, rStressVector.size() );
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // 2.-Compute e = 0.5*(inv(C) - I)
    const Matrix C_tensor = prod(trans(F),F);
    ConstitutiveLawUtilities<VoigtSize>::CalculateGreenLagrangianStrain(C_tensor, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicNeoHookean3D::CalculateAlmansiStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // 2.-COmpute e = 0.5*(1-inv(B))
    const Matrix B_tensor = prod(F,trans(F));
    ConstitutiveLawUtilities<VoigtSize>::CalculateAlmansiStrain(B_tensor, rStrainVector);
}

} // Namespace Kratos
