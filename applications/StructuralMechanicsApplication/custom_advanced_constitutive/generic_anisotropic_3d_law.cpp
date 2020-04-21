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
#include "generic_anisotropic_3d_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"


namespace Kratos
{
ConstitutiveLaw::Pointer GenericAnisotropic3DLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<GenericAnisotropic3DLaw>();
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_Cauchy);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_Kirchhoff);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // this strain in the real anisotropic space
    Vector& r_strain_vector = rValues.GetStrainVector();

    const Properties& r_material_properties = rValues.GetMaterialProperties();
   
    // We create the rValues for the isotropic CL
    const auto it_cl_begin                    = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl                = *(it_cl_begin);
    ConstitutiveLaw::Parameters values_iso_cl = rValues;
    values_iso_cl.SetMaterialProperties(r_props_iso_cl);

    // Here we rotate the orthotropy misalignment with respect to the global axes
    // {
    //}
    // todo ...

    // We compute the mappers As and Ae
    Matrix stress_mapper(VoigtSize, VoigtSize), strain_mapper(VoigtSize, VoigtSize);
    Matrix stress_mapper_inv(VoigtSize, VoigtSize); // The inverse of As
    Matrix isotropic_elastic_matrix(VoigtSize, VoigtSize), anisotropic_elastic_matrix(VoigtSize, VoigtSize);

    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rValues, stress_mapper, stress_mapper_inv);
    this->CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
    this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                isotropic_elastic_matrix, stress_mapper, 
                                                                                strain_mapper);
    // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal
    // TODO What it is F driven???????
    Vector &r_iso_strain_vector = values_iso_cl.GetStrainVector();
    r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector); // mapped
    mpIsotropicCL->CalculateMaterialResponsePK2(values_iso_cl);
    const Vector& r_iso_stress_vector = values_iso_cl.GetStressVector();

    // Finally, we map the stresses to the real space: Sreal = inv(As)Siso
    Vector &r_real_stress_vector = rValues.GetStressVector();
    r_real_stress_vector = prod(stress_mapper_inv, r_iso_stress_vector);

    // Here we revert the rotations
    // {
    //}
    // todo ...

} // End CalculateMaterialResponseCauchy


/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateElasticMatrix(
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

void GenericAnisotropic3DLaw::CalculateOrthotropicElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    if (rElasticityTensor.size1() != VoigtSize || rElasticityTensor.size2() != VoigtSize)
        rElasticityTensor.resize(VoigtSize, VoigtSize, false);
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
    rElasticityTensor(2, 2) = -E2 * E3 * P1 * (-E2 * v12 * v12 + E1);
    rElasticityTensor(3, 3) = (E2 * P2) / (P2 + v12 * (P2 + P3) + E1 * E2) / 2.0;
    rElasticityTensor(4, 4) = (E3 * P3) / (P3 + v23 * (P3 + P6) + E2 * E3) / 2.0;
    rElasticityTensor(5, 5) = (E3 * P2) / (P2 + v13 * (P2 + P6) + E1 * E3) / 2.0;
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // this strain in the real anisotropic space
    Vector& r_strain_vector = rValues.GetStrainVector();
    const Properties& r_material_properties = rValues.GetMaterialProperties();
   
    // We create the rValues for the isotropic CL
    const auto it_cl_begin                    = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl                = *(it_cl_begin);
    ConstitutiveLaw::Parameters values_iso_cl = rValues;
    values_iso_cl.SetMaterialProperties(r_props_iso_cl);

    // Here we rotate the orthotropy misalignment with respect to the global axes
    // {
    //}
    // todo ...

    // We compute the mappers As and Ae
    Matrix stress_mapper(VoigtSize, VoigtSize), strain_mapper(VoigtSize, VoigtSize);
    Matrix stress_mapper_inv(VoigtSize, VoigtSize); // The inverse of As
    Matrix isotropic_elastic_matrix(VoigtSize, VoigtSize), anisotropic_elastic_matrix(VoigtSize, VoigtSize);

    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rValues, stress_mapper, stress_mapper_inv);
    this->CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
    this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                isotropic_elastic_matrix, stress_mapper, 
                                                                                strain_mapper);
    // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal
    // TODO What it is F driven???????
    Vector &r_iso_strain_vector = values_iso_cl.GetStrainVector();
    r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector); // mapped
    mpIsotropicCL->FinalizeMaterialResponsePK2(values_iso_cl);
}

/***********************************************************************************/
/***********************************************************************************/

double& GenericAnisotropic3DLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& GenericAnisotropic3DLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& GenericAnisotropic3DLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropic3DLaw::Has(const Variable<bool>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropic3DLaw::Has(const Variable<double>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropic3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool GenericAnisotropic3DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

double& GenericAnisotropic3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& GenericAnisotropic3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        this->CalculateCauchyGreenStrain(rParameterValues, rValue);
        return rValue;
    } else if (rThisVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponsePK2(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;

    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;

    } else if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseKirchhoff(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;
        
    } else if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        if (mpIsotropicCL->Has(PLASTIC_STRAIN_VECTOR)) {

            const Vector &r_plastic_strain = mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
            const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        
            // We create the rParameterValues for the isotropic CL
            const auto it_cl_begin                    = r_material_properties.GetSubProperties().begin();
            const auto& r_props_iso_cl                = *(it_cl_begin);
            ConstitutiveLaw::Parameters values_iso_cl = rParameterValues;
            values_iso_cl.SetMaterialProperties(r_props_iso_cl);

            // We compute the mappers As and Ae
            Matrix stress_mapper(VoigtSize, VoigtSize), strain_mapper(VoigtSize, VoigtSize);
            Matrix stress_mapper_inv(VoigtSize, VoigtSize); // The inverse of As
            Matrix isotropic_elastic_matrix(VoigtSize, VoigtSize), anisotropic_elastic_matrix(VoigtSize, VoigtSize);

            ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rParameterValues, stress_mapper, stress_mapper_inv);
            this->CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
            this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
            ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                        isotropic_elastic_matrix, stress_mapper, 
                                                                                        strain_mapper);
            Matrix invAe(VoigtSize, VoigtSize);
            double aux_det;
            MathUtils<double>::InvertMatrix(strain_mapper, invAe, aux_det);
            rValue = prod(invAe, r_plastic_strain);
            return rValue;
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_isotropic_cl = *(it_cl_begin);
    KRATOS_ERROR_IF_NOT(r_props_isotropic_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    mpIsotropicCL = r_props_isotropic_cl[CONSTITUTIVE_LAW]->Clone();
    mpIsotropicCL->InitializeMaterial(r_props_isotropic_cl, rElementGeometry, rShapeFunctionsValues);

    // We check now the dimension of the CL pointer, must be 3D
    KRATOS_ERROR_IF_NOT(mpIsotropicCL->GetStrainSize() == 6) << "The slave CL has a dimension lower than 3, not possible" << std::endl;

    // Let's check variables
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_X))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_X not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Y))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Y not defined in properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XY)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XY not defined in properties" << std::endl;
    if (VoigtSize == 6) {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Z))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_Z not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XZ)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_XZ not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO_YZ)) << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO_YZ not defined in properties" << std::endl;
    }

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_X))  << "YOUNG_MODULUS_X not defined in properties"  << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_Y))  << "YOUNG_MODULUS_Y not defined in properties"  << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_XY)) << "POISSON_RATIO_XY not defined in properties" << std::endl;
    if (VoigtSize == 6) {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS_Z))  << "YOUNG_MODULUS_Z not defined in properties"  << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_XZ)) << "POISSON_RATIO_XZ not defined in properties" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO_YZ)) << "POISSON_RATIO_YZ not defined in properties" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& GenericAnisotropic3DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    mpIsotropicCL->InitializeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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

void GenericAnisotropic3DLaw::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    const SizeType space_dimension = this->WorkingSpaceDimension();

    // Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1()!= Dimension || F.size2() != Dimension)
        << "expected size of F " << Dimension << "x" << Dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

    Matrix E_tensor = prod(trans(F), F);
    for(unsigned int i = 0; i < Dimension; ++i)
        E_tensor(i, i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos
