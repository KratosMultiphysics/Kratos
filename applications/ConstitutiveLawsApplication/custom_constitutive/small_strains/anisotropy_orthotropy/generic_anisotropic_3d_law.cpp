// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
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

#include "constitutive_laws_application_variables.h"
#include "generic_anisotropic_3d_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
ConstitutiveLaw::Pointer GenericAnisotropicLaw<TDim>::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<GenericAnisotropicLaw>();
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector              = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f   = rValues.GetDeformationGradientF();
    const double determinant_f             = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_Cauchy);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f   = rValues.GetDeformationGradientF();
    const double determinant_f             = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_Kirchhoff);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateRotationMatrixVoigt(
    const Properties& rMaterialProperties,
    BoundedMatrixVoigtType &rT
)
{
    if (rMaterialProperties.Has(EULER_ANGLES) && MathUtils<double>::Norm3(rMaterialProperties[EULER_ANGLES]) > machine_tolerance) {
        BoundedMatrixType rotation_matrix;
        const Vector& r_euler_angles = rMaterialProperties[EULER_ANGLES];
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(r_euler_angles[0], r_euler_angles[1], r_euler_angles[2], rotation_matrix);
        ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(rotation_matrix, rT);
    } else {
        noalias(rT) = IdentityMatrix(VoigtSize, VoigtSize);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    if (!flag_strain) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        ConstitutiveLawUtilities<VoigtSize>::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // this strain in the real anisotropic space
    const Vector real_strain_vector         = rValues.GetStrainVector();
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // We create the rValues for the isotropic CL
    const auto it_cl_begin     = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl = *(it_cl_begin);
    rValues.SetMaterialProperties(r_props_iso_cl);

    if (flag_stress) {

        // Here we compute the rotation tensors due to the angles of the local and global axes
        BoundedMatrixVoigtType voigt_rotation_matrix;
        CalculateRotationMatrixVoigt(r_material_properties, voigt_rotation_matrix);

        // We compute the mappers As and Ae
        BoundedMatrixVoigtType stress_mapper, strain_mapper;
        BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
        BoundedMatrixVoigtType anisotropic_elastic_matrix;
        Matrix isotropic_elastic_matrix;

        this->CalculateAnisotropicStressMapperMatrix(r_material_properties, stress_mapper, stress_mapper_inv);
        mpIsotropicCL->CalculateValue(rValues, CONSTITUTIVE_MATRIX, isotropic_elastic_matrix); // Takes the props of the iso cl
        this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
        this->CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix, isotropic_elastic_matrix, stress_mapper, strain_mapper);

        Vector &r_iso_strain_vector = rValues.GetStrainVector();

        // Now we rotate the strain Eglob-> Eloc
        r_iso_strain_vector = prod((voigt_rotation_matrix), r_iso_strain_vector);

        // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal,loc
        r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector);

        // Integrate the isotropic constitutive law
        mpIsotropicCL->CalculateMaterialResponsePK2(rValues);
        Vector& r_iso_stress_vector = rValues.GetStressVector();

        // We map the stresses to the real space: Sreal = inv(As)Siso
        r_iso_stress_vector = prod(stress_mapper_inv, r_iso_stress_vector);

        // we return to the global coordinates the stress Sglob
        r_iso_stress_vector = prod(trans(voigt_rotation_matrix), r_iso_stress_vector);

        if (flag_const_tensor) {
            // Finally we map the tangent tensor: C_aniso = inv(As)*C_iso*Ae
            Matrix &r_anisotropic_tangent_matrix  = rValues.GetConstitutiveMatrix();
            r_anisotropic_tangent_matrix = prod(stress_mapper_inv, Matrix(prod(r_anisotropic_tangent_matrix, strain_mapper)));

            // Now we rotate to the global coordinates
            r_anisotropic_tangent_matrix = prod((trans(voigt_rotation_matrix)), Matrix(prod(r_anisotropic_tangent_matrix, (voigt_rotation_matrix))));
        }
        noalias(r_iso_strain_vector) = real_strain_vector;
    }
    rValues.SetMaterialProperties(r_material_properties);
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateOrthotropicElasticMatrix(
    BoundedMatrixVoigtType& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    if constexpr (Dimension == 3) {
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateOrthotropicElasticMatrix(rElasticityTensor, rMaterialProperties);
    } else {

    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    Flags& r_flags = rValues.GetOptions();
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        ConstitutiveLawUtilities<VoigtSize>::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // this strain in the real anisotropic space
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // We create the rValues for the isotropic CL
    const auto it_cl_begin     = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl = *(it_cl_begin);
    rValues.SetMaterialProperties(r_props_iso_cl);

    // Here we compute the rotation tensors due to the angles of the local and global axes
    BoundedMatrixVoigtType voigt_rotation_matrix;
    CalculateRotationMatrixVoigt(r_material_properties, voigt_rotation_matrix);

    // We compute the mappers As and Ae
    BoundedMatrixVoigtType stress_mapper, strain_mapper;
    BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
    BoundedMatrixVoigtType anisotropic_elastic_matrix;
    Matrix isotropic_elastic_matrix;

    this->CalculateAnisotropicStressMapperMatrix(r_material_properties, stress_mapper, stress_mapper_inv);
    mpIsotropicCL->CalculateValue(rValues, CONSTITUTIVE_MATRIX, isotropic_elastic_matrix); // takes the props of the iso cl
    this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
    this->CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix, isotropic_elastic_matrix, stress_mapper, strain_mapper);
    Vector &r_iso_strain_vector = rValues.GetStrainVector();
    // Now we rotate the strain Eglob-> Eloc
    r_iso_strain_vector = prod(voigt_rotation_matrix, r_iso_strain_vector);

    // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal,loc
    r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector); // mapped

    // Integrate the isotropic constitutive law
    mpIsotropicCL->FinalizeMaterialResponsePK2(rValues);

    rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double& GenericAnisotropicLaw<TDim>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
Vector& GenericAnisotropicLaw<TDim>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
Matrix& GenericAnisotropicLaw<TDim>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
bool GenericAnisotropicLaw<TDim>::Has(const Variable<bool>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
bool GenericAnisotropicLaw<TDim>::Has(const Variable<double>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
bool GenericAnisotropicLaw<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
bool GenericAnisotropicLaw<TDim>::Has(const Variable<Matrix>& rThisVariable)
{
    return mpIsotropicCL->Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double& GenericAnisotropicLaw<TDim>::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
Vector& GenericAnisotropicLaw<TDim>::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        Flags& r_flags = rParameterValues.GetOptions();
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            ConstitutiveLawUtilities<VoigtSize>::CalculateCauchyGreenStrain(rParameterValues, rValue);
        } else {
            noalias(rValue) = rParameterValues.GetStrainVector();
        }
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
        noalias(rValue) = rParameterValues.GetStressVector();

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
        noalias(rValue) = rParameterValues.GetStressVector();

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

            // Here we compute the rotation tensors due to the angles of the local and global axes
            BoundedMatrixVoigtType voigt_rotation_matrix, inv_voigt_rotation_matrix;
            CalculateRotationMatrixVoigt(r_material_properties, voigt_rotation_matrix);
            double det = 0.0;
            MathUtils<double>::InvertMatrix(voigt_rotation_matrix, inv_voigt_rotation_matrix, det);

            // We compute the mappers As and Ae
            BoundedMatrixVoigtType stress_mapper, strain_mapper;
            BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
            BoundedMatrixVoigtType anisotropic_elastic_matrix;
            Matrix isotropic_elastic_matrix;

            this->CalculateAnisotropicStressMapperMatrix(r_material_properties, stress_mapper, stress_mapper_inv);
            mpIsotropicCL->CalculateValue(values_iso_cl, CONSTITUTIVE_MATRIX, isotropic_elastic_matrix);
            this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
            this->CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix, isotropic_elastic_matrix, stress_mapper, strain_mapper);
            BoundedMatrixVoigtType invAe;
            double aux_det;
            MathUtils<double>::InvertMatrix(strain_mapper, invAe, aux_det);

            // We mapp to the local anisotropic space
            noalias(rValue) = prod(invAe, r_plastic_strain);

            // We rotate to the global system
            rValue = prod(inv_voigt_rotation_matrix, rValue);
            return rValue;
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateAnisotropicStressMapperMatrix(
    const Properties& rProperties,
    BoundedMatrixVoigtType& rAs,
    BoundedMatrixVoigtType& rAsInv
)
{
    rAs.clear();
    rAsInv.clear();
    const Vector &r_iso_aniso_yield_ratios = rProperties[ISOTROPIC_ANISOTROPIC_YIELD_RATIO];
    if constexpr (VoigtSize == 6) {
        rAs(0, 0) = r_iso_aniso_yield_ratios(0);
        rAs(1, 1) = r_iso_aniso_yield_ratios(1);
        rAs(2, 2) = r_iso_aniso_yield_ratios(2);
        rAs(3, 3) = r_iso_aniso_yield_ratios(3);
        rAs(4, 4) = r_iso_aniso_yield_ratios(4);
        rAs(5, 5) = r_iso_aniso_yield_ratios(5);
    } else {
        rAs(0, 0) = r_iso_aniso_yield_ratios(0);
        rAs(1, 1) = r_iso_aniso_yield_ratios(1);
        rAs(2, 2) = r_iso_aniso_yield_ratios(2);
    }

    for (IndexType i = 0; i < VoigtSize; ++i)
        rAsInv(i, i) = 1.0 / rAs(i, i);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateAnisotropicStrainMapperMatrix(
    const BoundedMatrixVoigtType& rAnisotropicElasticMatrix,
    const BoundedMatrixVoigtType& rIsotropicElasticMatrix,
    const BoundedMatrixVoigtType &rAs,
    BoundedMatrixVoigtType& rAe
)
{
    Matrix inv_isotropic_elastic_matrix(VoigtSize, VoigtSize);
    double aux_det;
    MathUtils<double>::InvertMatrix(rIsotropicElasticMatrix, inv_isotropic_elastic_matrix, aux_det);
    noalias(rAe) = prod(inv_isotropic_elastic_matrix, Matrix(prod(rAs, rAnisotropicElasticMatrix)));
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_isotropic_cl = *(it_cl_begin);
    KRATOS_ERROR_IF_NOT(r_props_isotropic_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    mpIsotropicCL = r_props_isotropic_cl[CONSTITUTIVE_LAW]->Clone();
    mpIsotropicCL->InitializeMaterial(r_props_isotropic_cl, rElementGeometry, rShapeFunctionsValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
Matrix& GenericAnisotropicLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return mpIsotropicCL->CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::InitializeMaterialResponsePK2(Parameters& rValues)
{
    mpIsotropicCL->InitializeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
int GenericAnisotropicLaw<TDim>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto it_cl_begin     = rMaterialProperties.GetSubProperties().begin();
    const auto& r_props_iso_cl = *(it_cl_begin);

    // We check now the dimension of the CL pointer, must be 3D
    KRATOS_ERROR_IF_NOT(mpIsotropicCL->GetStrainSize() == VoigtSize) << "The slave CL has a different dimension of the Generic Anisotropic CL..." << std::endl;

    // Let's check variables
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO not defined in properties" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ORTHOTROPIC_ELASTIC_CONSTANTS)) << "The ORTHOTROPIC_ELASTIC_CONSTANTS are not defined" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[ORTHOTROPIC_ELASTIC_CONSTANTS].size() == VoigtSize) << "The dimension of the ORTHOTROPIC_ELASTIC_CONSTANTS is incorrect" << std::endl;


    return mpIsotropicCL->Check(r_props_iso_cl, rElementGeometry, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenericAnisotropicLaw<TDim>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_PK2, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_PK2, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericAnisotropicLaw<2>;
template class GenericAnisotropicLaw<3>;

} // namespace Kratos
