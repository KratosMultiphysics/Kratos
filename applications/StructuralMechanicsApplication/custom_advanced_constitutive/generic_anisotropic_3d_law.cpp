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

    Vector& r_stress_vector              = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, deformation_gradient_f, determinant_f,
                            ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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

void GenericAnisotropic3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
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

void GenericAnisotropic3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // this strain in the real anisotropic space
    const Vector real_strain_vector         = rValues.GetStrainVector();
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    if (!flag_strain) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        this->CalculateCauchyGreenStrain(rValues, r_strain_vector);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // We create the rValues for the isotropic CL
    const auto it_cl_begin                    = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl                = *(it_cl_begin);
    ConstitutiveLaw::Parameters values_iso_cl = rValues;
    values_iso_cl.SetMaterialProperties(r_props_iso_cl);

    if (flag_stress) {

        // Here we compute the rotation tensors due to the angles of the local and global axes
        BoundedMatrixType rotation_matrix;
        BoundedMatrixVoigtType voigt_rotation_matrix, inv_voigt_rotation_matrix;

        if (r_material_properties.Has(EULER_ANGLE_PHI)   &&
            r_material_properties.Has(EULER_ANGLE_THETA) &&
            r_material_properties.Has(EULER_ANGLE_HI)    &&
            std::abs(r_material_properties[EULER_ANGLE_PHI]) + 
            std::abs(r_material_properties[EULER_ANGLE_THETA]) + 
            std::abs(r_material_properties[EULER_ANGLE_HI]) > machine_tolerance) {
                ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(
                    r_material_properties[EULER_ANGLE_PHI], r_material_properties[EULER_ANGLE_THETA],
                    r_material_properties[EULER_ANGLE_HI], rotation_matrix);
                ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(
                    (rotation_matrix),
                    voigt_rotation_matrix);
            double aux_det = 0.0;
            MathUtils<double>::InvertMatrix(voigt_rotation_matrix, inv_voigt_rotation_matrix, aux_det);
        } else {
            noalias(rotation_matrix)           = IdentityMatrix(Dimension, Dimension);
            noalias(voigt_rotation_matrix)     = IdentityMatrix(VoigtSize, VoigtSize);
            noalias(inv_voigt_rotation_matrix) = IdentityMatrix(VoigtSize, VoigtSize);
        }
        
        // We compute the mappers As and Ae
        BoundedMatrixVoigtType stress_mapper, strain_mapper;
        BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
        BoundedMatrixVoigtType isotropic_elastic_matrix, anisotropic_elastic_matrix;

        ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rValues, stress_mapper, stress_mapper_inv);
        ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
        this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
        ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                    isotropic_elastic_matrix, stress_mapper, 
                                                                                    strain_mapper);
        Vector &r_iso_strain_vector = values_iso_cl.GetStrainVector();
        // Now we rotate the strain Eglob-> Eloc
        r_iso_strain_vector = prod((voigt_rotation_matrix), r_iso_strain_vector);

        // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal,loc
        r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector);

        // Integrate the isotropic constitutive law
        mpIsotropicCL->CalculateMaterialResponsePK2(values_iso_cl);
        const Vector& r_iso_stress_vector = values_iso_cl.GetStressVector();

        // We map the stresses to the real space: Sreal = inv(As)Siso
        Vector &r_real_stress_vector = rValues.GetStressVector();
        noalias(r_real_stress_vector) = prod(stress_mapper_inv, r_iso_stress_vector);

        // we return to the global coordinates the stress Sglob
        r_real_stress_vector = prod(trans(voigt_rotation_matrix), r_real_stress_vector);

        if (flag_const_tensor) {
            // Finally we map the tangent tensor: C_aniso = inv(As)*C_iso*Ae
            Matrix &r_anisotropic_tangent_matrix  = rValues.GetConstitutiveMatrix();
            const Matrix& r_isotropic_tangent     = values_iso_cl.GetConstitutiveMatrix();
            noalias(r_anisotropic_tangent_matrix) = prod(stress_mapper_inv, Matrix(prod(r_isotropic_tangent, strain_mapper)));

            // Now we rotate to the global coordinates
            noalias(r_anisotropic_tangent_matrix) = prod((trans(voigt_rotation_matrix)), Matrix(prod(r_isotropic_tangent, (voigt_rotation_matrix))));
        }  
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateOrthotropicElasticMatrix(
    BoundedMatrixVoigtType& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    KRATOS_TRY

    if (rElasticityTensor.size1() != VoigtSize || rElasticityTensor.size2() != VoigtSize)
        rElasticityTensor.resize(VoigtSize, VoigtSize, false);
    noalias(rElasticityTensor) = ZeroMatrix(VoigtSize, VoigtSize);

    const double Ex  = rMaterialProperties[YOUNG_MODULUS_X];
    const double Ey  = rMaterialProperties[YOUNG_MODULUS_Y];
    const double Ez  = rMaterialProperties[YOUNG_MODULUS_Z];
    const double vxy = rMaterialProperties[POISSON_RATIO_XY];
    const double vyz = rMaterialProperties[POISSON_RATIO_YZ];
    const double vxz = rMaterialProperties[POISSON_RATIO_XZ];

    const double vyx = vxy * Ey / Ex;
    const double vzx = vxz * Ez / Ex;
    const double vzy = vyz * Ez / Ey;

    KRATOS_ERROR_IF(vyx > 0.5) << "The Poisson_yx is greater than 0.5." << std::endl;
    KRATOS_ERROR_IF(vzx > 0.5) << "The Poisson_zx is greater than 0.5." << std::endl;
    KRATOS_ERROR_IF(vzy > 0.5) << "The Poisson_zy is greater than 0.5." << std::endl;

    const double Gxy   = 1.0 / ((1.0 + vyx) / Ex + (1.0 + vxy) / Ey);
    const double Gxz   = 1.0 / ((1.0 + vzx) / Ex + (1.0 + vxz) / Ez);
    const double Gyz   = 1.0 / ((1.0 + vzy) / Ey + (1.0 + vyz) / Ez);

    const double ctant = 1.0 / (1.0 - vxy * vyx - vzy * vyz - vzx * vxz - vxy * vyz * vzx - vxz * vyx * vzy);

    rElasticityTensor(0, 0) = Ex * (1.0 - vyz * vzy) * ctant;
    rElasticityTensor(0, 1) = Ex * (vyx + vyz * vzx) * ctant;
    rElasticityTensor(1, 0) = Ey * (vxy + vzy * vxz) * ctant;

    rElasticityTensor(0, 2) = Ex * (vxz + vyx * vzy) * ctant;
    rElasticityTensor(2, 0) = Ez * (vxz + vxy * vyz) * ctant;
    rElasticityTensor(1, 1) = Ey * (1.0 - vxz * vzx) * ctant;

    rElasticityTensor(1, 2) = Ey * (vzy + vxy * vzx) * ctant;
    rElasticityTensor(2, 1) = Ez * (vyz + vyx * vxz) * ctant;
    rElasticityTensor(2, 2) = Ez * (1.0 - vxy * vyx) * ctant;

    rElasticityTensor(3, 3) = Gxy;
    rElasticityTensor(4, 4) = Gyz;
    rElasticityTensor(5, 5) = Gxz;

    KRATOS_CATCH("")
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
    const Vector real_strain_vector         = rValues.GetStrainVector();
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // We create the rValues for the isotropic CL
    const auto it_cl_begin                    = r_material_properties.GetSubProperties().begin();
    const auto& r_props_iso_cl                = *(it_cl_begin);
    ConstitutiveLaw::Parameters values_iso_cl = rValues;
    values_iso_cl.SetMaterialProperties(r_props_iso_cl);

    Flags& r_flags = rParameterValues.GetOptions();
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        this->CalculateCauchyGreenStrain(rValues, r_strain_vector);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // Here we compute the rotation tensors due to the angles of the local and global axes
    BoundedMatrixType rotation_matrix;
    BoundedMatrixVoigtType voigt_rotation_matrix;

    if (r_material_properties.Has(EULER_ANGLE_PHI)   &&
        r_material_properties.Has(EULER_ANGLE_THETA) &&
        r_material_properties.Has(EULER_ANGLE_HI)    &&
        std::abs(r_material_properties[EULER_ANGLE_PHI]) + 
        std::abs(r_material_properties[EULER_ANGLE_THETA]) + 
        std::abs(r_material_properties[EULER_ANGLE_HI]) > machine_tolerance) {
            ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(
                r_material_properties[EULER_ANGLE_PHI], r_material_properties[EULER_ANGLE_THETA],
                r_material_properties[EULER_ANGLE_HI], rotation_matrix);
            ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(
                (rotation_matrix),
                voigt_rotation_matrix);
    } else {
        noalias(rotation_matrix)           = IdentityMatrix(Dimension, Dimension);
        noalias(voigt_rotation_matrix)     = IdentityMatrix(VoigtSize, VoigtSize);
    }
    
    // We compute the mappers As and Ae
    BoundedMatrixVoigtType stress_mapper, strain_mapper;
    BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
    BoundedMatrixVoigtType isotropic_elastic_matrix, anisotropic_elastic_matrix;

    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rValues, stress_mapper, stress_mapper_inv);
    ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
    this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
    ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                isotropic_elastic_matrix, stress_mapper, 
                                                                                strain_mapper);
    Vector &r_iso_strain_vector = values_iso_cl.GetStrainVector();
    // Now we rotate the strain Eglob-> Eloc
    r_iso_strain_vector = prod((voigt_rotation_matrix), r_iso_strain_vector);

    // Now we map the strains to the isotropic fictitious space: Eiso = Ae*Ereal,loc
    r_iso_strain_vector = prod(strain_mapper, r_iso_strain_vector); // mapped

    // Integrate the isotropic constitutive law
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
    return false;
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
        Flags& r_flags = rParameterValues.GetOptions();
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateCauchyGreenStrain(rParameterValues, rValue);
        } else {
            rValue = rParameterValues.GetStrainVector();
        }
        KRATOS_WATCH(rValue)
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
        KRATOS_WATCH(rValue)
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
            BoundedMatrixType rotation_matrix;
            BoundedMatrixVoigtType voigt_rotation_matrix, inv_voigt_rotation_matrix;

            if (r_material_properties.Has(EULER_ANGLE_PHI)   &&
                r_material_properties.Has(EULER_ANGLE_THETA) &&
                r_material_properties.Has(EULER_ANGLE_HI)    &&
                std::abs(r_material_properties[EULER_ANGLE_PHI]) + 
                std::abs(r_material_properties[EULER_ANGLE_THETA]) + 
                std::abs(r_material_properties[EULER_ANGLE_HI]) > machine_tolerance) {
                    ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(
                        r_material_properties[EULER_ANGLE_PHI], r_material_properties[EULER_ANGLE_THETA],
                        r_material_properties[EULER_ANGLE_HI], rotation_matrix);
                    ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(
                        (rotation_matrix),
                        voigt_rotation_matrix);
                double aux_det = 0.0;
                MathUtils<double>::InvertMatrix(voigt_rotation_matrix, inv_voigt_rotation_matrix, aux_det);
            } else {
                noalias(rotation_matrix)           = IdentityMatrix(Dimension, Dimension);
                noalias(voigt_rotation_matrix)     = IdentityMatrix(VoigtSize, VoigtSize);
                noalias(inv_voigt_rotation_matrix) = IdentityMatrix(VoigtSize, VoigtSize);
            }

            // We compute the mappers As and Ae
            BoundedMatrixVoigtType stress_mapper, strain_mapper;
            BoundedMatrixVoigtType stress_mapper_inv; // The inverse of As
            BoundedMatrixVoigtType isotropic_elastic_matrix, anisotropic_elastic_matrix;

            ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStressMapperMatrix(rParameterValues, stress_mapper, stress_mapper_inv);
            ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(isotropic_elastic_matrix, r_props_iso_cl); // takes the props of the iso cl
            this->CalculateOrthotropicElasticMatrix(anisotropic_elastic_matrix, r_material_properties);
            ConstitutiveLawUtilities<VoigtSize>::CalculateAnisotropicStrainMapperMatrix(anisotropic_elastic_matrix,
                                                                                        isotropic_elastic_matrix, stress_mapper, 
                                                                                        strain_mapper);
            BoundedMatrixVoigtType invAe;
            double aux_det;
            MathUtils<double>::InvertMatrix(strain_mapper, invAe, aux_det);

            // We mapp to the local anisotropic space
            rValue = prod(invAe, r_plastic_strain);

            // We rotate to the global system
            rValue = prod(inv_voigt_rotation_matrix, r_plastic_strain);
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
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(ISOTROPIC_ANISOTROPIC_YIELD_RATIO))  << "ISOTROPIC_ANISOTROPIC_YIELD_RATIO not defined in properties" << std::endl;
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

void GenericAnisotropic3DLaw::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    // Compute total deformation gradient
    const BoundedMatrixType& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1()!= Dimension || F.size2() != Dimension)
        << "expected size of F " << Dimension << "x" << Dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

    BoundedMatrixType E_tensor = prod(trans(F), F);
    for(unsigned int i = 0; i < Dimension; ++i)
        E_tensor(i, i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

/***********************************************************************************/
/***********************************************************************************/

int GenericAnisotropic3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto it_cl_begin     = rMaterialProperties.GetSubProperties().begin();
    const auto& r_props_iso_cl = *(it_cl_begin);
    return mpIsotropicCL->Check(r_props_iso_cl, rElementGeometry, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropic3DLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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

} // namespace Kratos
