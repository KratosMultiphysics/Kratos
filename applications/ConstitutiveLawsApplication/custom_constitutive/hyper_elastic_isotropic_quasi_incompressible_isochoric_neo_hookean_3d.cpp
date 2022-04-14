// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/hyper_elastic_isotropic_quasi_incompressible_isochoric_neo_hookean_3d.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D()
    : HyperElasticIsotropicNeoHookean3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D(const HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D& rOther)
    : HyperElasticIsotropicNeoHookean3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::~HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f             = rValues.GetDeterminantF();

    TransformStresses(r_stress_vector, r_deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& r_strain_vector                  = rValues.GetStrainVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    Matrix C_tensor(dimension, dimension), inverse_C_tensor(dimension, dimension);

    if(r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        BaseType::CalculateGreenLagrangianStrain(rValues, r_strain_vector);
        noalias(C_tensor) = prod(trans(r_deformation_gradient_f), r_deformation_gradient_f);
    } else {
        Matrix strain_tensor(dimension, dimension);
        noalias(strain_tensor) = MathUtils<double>::StrainVectorToTensor(r_strain_vector);
        noalias(C_tensor) = 2.0 * strain_tensor + IdentityMatrix(dimension);
        determinant_f = std::sqrt(MathUtils<double>::Det(C_tensor));
    }

    double aux_det;
    MathUtils<double>::InvertMatrix(C_tensor, inverse_C_tensor, aux_det);

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{


}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);

    Vector& r_stress_vector       = rValues.GetStressVector();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f    = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    r_stress_vector       /= determinant_f;
    r_constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = VoigtSize;

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

int HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
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

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateConstitutiveMatrixPK2(
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

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateConstitutiveMatrixKirchhoff(
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

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculatePK2Stress(
    const Matrix& rInvCTensor,
    Vector& rStressVector,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    const SizeType dimension = WorkingSpaceDimension();
    Matrix stress_matrix(dimension, dimension);
    const Matrix Id = IdentityMatrix(dimension);
    noalias(stress_matrix) = LameLambda * std::log(DeterminantF) * rInvCTensor + LameMu * ( Id - rInvCTensor );
    noalias(rStressVector) = MathUtils<double>::StressTensorToVector( stress_matrix, GetStrainSize());
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::CalculateKirchhoffStress(
    const Matrix& rBTensor,
    Vector& rStressVector,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu
    )
{
    const SizeType dimension = WorkingSpaceDimension();
    Matrix stress_matrix(dimension, dimension);
    const Matrix Id = IdentityMatrix(dimension);
    noalias(stress_matrix) = LameLambda * std::log(DeterminantF) * Id + LameMu * (rBTensor - Id);
    noalias(rStressVector) = MathUtils<double>::StressTensorToVector(stress_matrix, GetStrainSize());
}

} // Namespace Kratos
