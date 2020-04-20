// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "generic_anisotropic_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"


namespace Kratos
{
ConstitutiveLaw::Pointer GenericAnisotropicLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<GenericAnisotropicLaw>();
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // // Some auxiliar values
    // const SizeType dimension = WorkingSpaceDimension();
    // const SizeType voigt_size = GetStrainSize();

    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();

    // // Previous flags saved
    // const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    // const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    // const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // // The deformation gradient
    // if (rValues.IsSetDeterminantF()) {
    //     const double determinant_f = rValues.GetDeterminantF();
    //     KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    // }
    // // In case the element has not computed the Strain
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     Vector& r_strain_vector = rValues.GetStrainVector();

    //     Matrix F_deformation_gradient(dimension, dimension);
    //     this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
    //     const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
    //     // Doing resize in case is needed
    //     if (r_strain_vector.size() != voigt_size)
    //         r_strain_vector.resize(voigt_size);

    //      // Identity matrix
    //     Matrix identity_matrix(dimension, dimension);
    //     for (IndexType i = 0; i < dimension; ++i) {
    //         for (IndexType j = 0; j < dimension; ++j) {
    //             if (i == j) identity_matrix(i, j) = 1.0;
    //             else identity_matrix(i, j) = 0.0;
    //         }
    //     }

    //     // Calculating the inverse of the left Cauchy tensor
    //     Matrix inverse_B_tensor(dimension, dimension);
    //     double aux_det_b = 0;
    //     MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

    //     // Calculate E matrix
    //     const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
    //     // Almansi Strain Calculation
    //     r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    // }

    // if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
    //     // Set new flags
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    //     // Total strain vector
    //     Vector& r_strain_vector = rValues.GetStrainVector();
    //     Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
    //     Vector fiber_stress_vector, matrix_stress_vector;
    //     this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
    //                                                 fiber_stress_vector,
    //                                                 matrix_stress_vector,
    //                                                 r_material_properties,
    //                                                 rValues,
    //                                                 serial_strain_matrix_old);
    //     Vector& r_integrated_stress_vector = rValues.GetStressVector();
    //     noalias(r_integrated_stress_vector) = mFiberVolumetricParticipation * fiber_stress_vector 
    //                                  + (1.0 - mFiberVolumetricParticipation) * matrix_stress_vector;

    //     // Previous flags restored
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

    //     if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //         this->CalculateTangentTensor(rValues);
    //     }
    // }

} // End CalculateMaterialResponseCauchy


/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const double mu = E / (2.0 + 2.0 * poisson_ratio);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2.0 * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2.0 * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2.0 * mu;
    rElasticityTensor(3, 3) = mu;
    rElasticityTensor(4, 4) = mu;
    rElasticityTensor(5, 5) = mu;
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // const Vector& r_strain_vector = rValues.GetStrainVector();
    // mPreviousStrainVector = r_strain_vector;

    // // Recalculation to obtain the serial_strain_matrix and store the value
    // const SizeType voigt_size = GetStrainSize();

    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();

    // // Previous flags saved
    // const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    // const Properties& r_material_properties = rValues.GetMaterialProperties();

    // if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
    //     // Set new flags
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    //     // Total strain vector
    //     Vector& r_strain_vector = rValues.GetStrainVector();
    //     Vector fiber_stress_vector, matrix_stress_vector;
    //     this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
    //                                                 fiber_stress_vector,
    //                                                 matrix_stress_vector,
    //                                                 r_material_properties,
    //                                                 rValues,
    //                                                 mPreviousSerialStrainMatrix);
    //     // Previous flags restored
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    //     r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

    //     // We call the FinalizeMaterialResponse of the matrix and fiber CL
    //     auto& r_material_properties = rValues.GetMaterialProperties();
    //     const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    //     const auto& r_props_matrix_cl = *(it_cl_begin);
    //     const auto& r_props_fiber_cl = *(it_cl_begin + 1);
        
    //     ConstitutiveLaw::Parameters values_fiber  = rValues;
    //     ConstitutiveLaw::Parameters values_matrix = rValues;

    //     values_matrix.SetMaterialProperties(r_props_matrix_cl);
    //     values_fiber.SetMaterialProperties(r_props_fiber_cl);

    //     Matrix parallel_projector, serial_projector;
    //     this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
    //     Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

    //     this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, 
    //                                           mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector);

    //     values_fiber.SetStrainVector(fiber_strain_vector);
    //     values_matrix.SetStrainVector(matrix_strain_vector);

    //     mpMatrixConstitutiveLaw->FinalizeMaterialResponseCauchy(values_matrix);
    //     mpFiberConstitutiveLaw ->FinalizeMaterialResponseCauchy(values_fiber);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

double& GenericAnisotropicLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& GenericAnisotropicLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& GenericAnisotropicLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropicLaw::Has(const Variable<bool>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropicLaw::Has(const Variable<double>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropicLaw::Has(const Variable<Vector>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropicLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

double& GenericAnisotropicLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& GenericAnisotropicLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_isotropic_cl = *(it_cl_begin);
    KRATOS_ERROR_IF_NOT(r_props_isotropic_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    mpIsotropicCL = r_props_isotropic_cl[CONSTITUTIVE_LAW]->Clone();
    mpIsotropicCL->InitializeMaterial(r_props_isotropic_cl, rElementGeometry, rShapeFunctionsValues);

    // Let's check variables
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_X))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_X not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Y))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Y not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XY)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XY not defined in properties" << std::endl;
    if (mpIsotropicCL->GetStrainSize() == 6) {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Z))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Z not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XZ)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XZ not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_YZ)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_YZ not defined in properties" << std::endl;
    }

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_X))  << "YOUNG_MODULUS_X not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_Y))  << "YOUNG_MODULUS_Y not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_XY)) << "POISSON_RATIO_XY not defined in properties" << std::endl;
    if (mpIsotropicCL->GetStrainSize() == 6) {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_Z))  << "YOUNG_MODULUS_Z not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_XZ)) << "POISSON_RATIO_XZ not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_YZ)) << "POISSON_RATIO_YZ not defined in properties" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& GenericAnisotropicLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    mpIsotropicCL->InitializeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}
/***********************************************************************************/
/***********************************************************************************/
} // namespace Kratos
