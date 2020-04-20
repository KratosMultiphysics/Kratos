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

void GenericAnisotropicLaw::CalculateOrthotropicElasticMatrix( // TODO generalize for 2D
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

	const double E1  = rMaterialProperties[YOUNG_MODULUS_X];
	const double E2  = rMaterialProperties[YOUNG_MODULUS_Y];
	const double E3  = rMaterialProperties[YOUNG_MODULUS_Z];
	const double v12 = rMaterialProperties[POISSON_RATIO_XY];
	const double v23 = rMaterialProperties[POISSON_RATIO_YZ];
	const double v13 = rMaterialProperties[POISSON_RATIO_XZ];
    const double P1  = 1.0 / (E2 * E2 * v12 * v12 + 2.0 * E3 * E2 * v12 * v13 * v23 + E3 * E2 * v13 * v13 - E1 * E2 + E1 * E3 * v23 * v23);
    const double P2  = E1 * E1;
    const double P3  = E2 * E2;
    const double P4  = E1 * v23 + E2 * v12 * v13;
    const double P5  = E2 * v12 + E3 * v13 * v23;
    const double P6  = E3 * E3;

    rElasticityTensor(0, 0) = -P1 * P2 * (-E3 * v23 * v23 + E2);
    rElasticityTensor(0, 1) = -E1 * E2 * P1 * P5;
    rElasticityTensor(0, 2) = -E2 * E3 * P1 * (E1 * v13 + E1 * v12 * v23);
    rElasticityTensor(1, 0) = -E1 * E2 * P1 * P5;
    rElasticityTensor(1, 1) = -P1 * P3 * (-E3 * v13 * v13 + E1);
    rElasticityTensor(1, 2) = -E2 * E3 * P1 * P4;
    rElasticityTensor(2, 0) = -E1 * E2 * E3 * P1 * (v13 + v12 * v23);
    rElasticityTensor(2, 1) = -E2 * E3 * P1 * P4;
    rElasticityTensor(2, 2) = -E2*E3*P1*(- E2*v12*v12 + E1);
    rElasticityTensor(3, 3) = (E2 * P2) / (P2 + v12 * (P2 + P3) + E1 * E2) / 2.0;
    rElasticityTensor(4, 4) = (E3 * P3) / (P3 + v23 * (P3 + P6) + E2 * E3) / 2.0;
    rElasticityTensor(5, 5) = (E3 * P2) / (P2 + v13 * (P2 + P6) + E1 * E3) / 2.0;
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
    Vector& r_strain_vector = rValues.GetStrainVector();
    
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    


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
    // const Properties& r_material_properties = rValues.GetMaterialProperties();

    // const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    // const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    // if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
    //     KRATOS_ERROR << "Analytic solution not available" << std::endl;
    // } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
    //     // Calculates the Tangent Constitutive Tensor by perturbation (first order)
    //     TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    // } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
    //     // Calculates the Tangent Constitutive Tensor by perturbation (second order)
    //     TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    // }
}
/***********************************************************************************/
/***********************************************************************************/
} // namespace Kratos
