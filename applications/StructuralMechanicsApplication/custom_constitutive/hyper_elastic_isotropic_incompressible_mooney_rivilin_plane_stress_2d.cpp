// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/hyper_elastic_isotropic_incompressible_mooney_rivilin_plane_stress_2d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D(const HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::~HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D(){};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    CalculateMaterialResponsePK2(rValues);

    Vector &stress_vector = rValues.GetStressVector();
    const Matrix &deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();

    TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags &r_flags = rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();

    const Properties &material_properties = rValues.GetMaterialProperties();
    Vector &strain_vector = rValues.GetStrainVector();

    // The material properties
    const double c10 = material_properties[YOUNG_MODULUS]; // nav have to change the variable name
    const double c01 = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix &deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // We compute the right Cauchy-Green tensor (C):
    Matrix C_tensor = prod(trans(deformation_gradient_f), deformation_gradient_f);
    Matrix C_tensor_reduced(dimension, dimension);
    for (IndexType i = 0; i < dimension; ++i)
    {
        for (IndexType j = 0; j < dimension; ++j)
        {
            C_tensor_reduced(i, j) = C_tensor(i, j);
        }
    }

    // Inverse of the right Cauchy-Green tensor (C):
    double det_C_reduced;
    Matrix inverse_C_tensor_reduced(dimension, dimension);
    MathUtils<double>::InvertMatrix(C_tensor_reduced, inverse_C_tensor_reduced, det_C_reduced);

    double first_invariant_reduced = 0.0;

    for (IndexType i = 0; i < dimension; ++i)
    {
        first_invariant_reduced += C_tensor_reduced(i, i);
    }

    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
    {

        this->CalculateGreenLagrangianStrain(rValues, strain_vector);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {

        Matrix &constitutive_matrix = rValues.GetConstitutiveMatrix();

        CalculateConstitutiveMatrixPK2(constitutive_matrix, inverse_C_tensor_reduced, det_C_reduced, first_invariant_reduced, c10, c01);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {

        Vector &stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(C_tensor_reduced, inverse_C_tensor_reduced, stress_vector, det_C_reduced, first_invariant_reduced, c10, c01);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    // Get Values to compute the constitutive law:
    Flags &r_flags = rValues.GetOptions();

    const Properties &material_properties = rValues.GetMaterialProperties();
    Vector &strain_vector = rValues.GetStrainVector();
    Vector &stress_vector = rValues.GetStressVector();

    // The material properties
    const double c10 = material_properties[YOUNG_MODULUS]; // nav have to change the variable name
    const double c01 = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix &deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    if (r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
    {
        CalculateAlmansiStrain(rValues, strain_vector);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        Matrix &constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrixKirchhoff(constitutive_matrix, determinant_f, c10, c01);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        // We compute the left Cauchy-Green tensor (B):
        const Matrix B_tensor = prod(deformation_gradient_f, trans(deformation_gradient_f));
        CalculateKirchhoffStress(B_tensor, stress_vector, determinant_f, c10, c01);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);

    Vector &stress_vector = rValues.GetStressVector();
    Matrix &constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    stress_vector /= determinant_f;
    constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK1(rValues);
    //     rValues.Reset(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rValues);
    //     rValues.Reset(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseCauchy(rValues);
    //     rValues.Reset(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseKirchhoff(rValues);
    //     rValues.Reset(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK1(rValues);
    //     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rValues);
    //     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseCauchy(rValues);
    //     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    //     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    //     HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseKirchhoff(rValues);
    //     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

double &HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<double> &rThisVariable,
    double &rValue)
{
    const Properties &material_properties = rParameterValues.GetMaterialProperties();

    // The LAME parameters
    const double c10 = material_properties[YOUNG_MODULUS]; // nav have to change the variable name
    const double c01 = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix &deformation_gradient_f = rParameterValues.GetDeformationGradientF();

    // We compute the right Cauchy-Green tensor (C):
    const Matrix C_tensor = prod(trans(deformation_gradient_f), deformation_gradient_f);

    if (rThisVariable == STRAIN_ENERGY)
    {

        double first_invariant = 0.0;
        double second_invariant = 0.0;

        for (IndexType i = 0; i < C_tensor.size1(); ++i)
        {
            first_invariant += C_tensor(i, i);
        }

        for (IndexType i = 0; i < C_tensor.size1(); ++i)
        {
            for (IndexType j = 0; j < C_tensor.size1(); ++j)
            {
                second_invariant -= C_tensor(i, j) * C_tensor(i, j);
            }
        }

        second_invariant += first_invariant * first_invariant;
        second_invariant /= 2;
        rValue = c10 * (first_invariant - 3.0) + c01 * (first_invariant - 3.0);
    }

    return (rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector &HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<Vector> &rThisVariable,
    Vector &rValue)
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == HENCKY_STRAIN_VECTOR ||
        rThisVariable == BIOT_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR)
    {

        // Get Values to compute the constitutive law:
        Flags &r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

        // We compute the strain
        if (rThisVariable == STRAIN)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        }
        else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rParameterValues);
        }
        else if (rThisVariable == ALMANSI_STRAIN_VECTOR)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseKirchhoff(rParameterValues);
        }
        else if (rThisVariable == HENCKY_STRAIN_VECTOR)
        {
            const Matrix &deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans(deformation_gradient_f), deformation_gradient_f);
            Vector &r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_tensor, r_strain_vector);
        }
        else if (rThisVariable == BIOT_STRAIN_VECTOR)
        {
            const Matrix &deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans(deformation_gradient_f), deformation_gradient_f);
            Vector &r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateBiotStrain(C_tensor, r_strain_vector);
        }

        rValue = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }
    else if (rThisVariable == STRESSES ||
             rThisVariable == CAUCHY_STRESS_VECTOR ||
             rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
             rThisVariable == PK2_STRESS_VECTOR)
    {

        // Get Values to compute the constitutive law:
        Flags &r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        if (rThisVariable == STRESSES)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        }
        if (rThisVariable == KIRCHHOFF_STRESS_VECTOR)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseKirchhoff(rParameterValues);
        }
        if (rThisVariable == CAUCHY_STRESS_VECTOR)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponseCauchy(rParameterValues);
        }
        if (rThisVariable == PK2_STRESS_VECTOR)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }

    return (rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix &HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<Matrix> &rThisVariable,
    Matrix &rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF)
    {
        // Get Values to compute the constitutive law:
        Flags &r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

        // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        }
        else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rParameterValues);
        }
        else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF)
        {
            HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }

    return (rValue);
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::GetLawFeatures(Features &rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set(PLANE_STRESS_LAW);
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

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

int HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::Check(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "C10 is invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] <= -1e-9) << "C01 is invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is invalid value " << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateConstitutiveMatrixPK2(
    Matrix &rConstitutiveMatrix,
    const Matrix &rInverseCTensor,
    const double DeterminantC,
    const double FirstInvariant,
    const double C10,
    const double C01)
{
    rConstitutiveMatrix.clear();
    double InvCdotInvC, InvCtimesInvC, ItimesI, DC_DC, ItimesInvC, InvCtimesI;
    const SizeType dimension = WorkingSpaceDimension();
    Matrix I = IdentityMatrix(dimension, dimension);

    for (IndexType i = 0; i < 3; ++i)
    {
        const IndexType i0 = this->msIndexVoigt2D3C[i][0];
        const IndexType i1 = this->msIndexVoigt2D3C[i][1];

        for (IndexType j = 0; j < 3; ++j)
        {
            const IndexType j0 = this->msIndexVoigt2D3C[j][0];
            const IndexType j1 = this->msIndexVoigt2D3C[j][1];

            InvCdotInvC = 0.5 * (rInverseCTensor(i0, j0) * rInverseCTensor(i1, j1) + rInverseCTensor(i0, j1) * rInverseCTensor(i1, j0)); // sym
            //InvCdotInvC = rInverseCTensor(i0, j0) * rInverseCTensor(i1, j1);
            InvCtimesInvC = rInverseCTensor(i0, i1) * rInverseCTensor(j0, j1);
            ItimesI = I(i0, i1) * I(j0, j1);
            ItimesInvC = I(i0, i1) * rInverseCTensor(j0, j1);
            InvCtimesI = rInverseCTensor(i0, i1) * I(j0, j1);
            DC_DC = 0.5 * (I(i0, j0) * I(i1, j1) + I(i0, j1) * I(i1, j0)); // sym
            //DC_DC = I(i0, j0) * I(i1, j1) ;

              rConstitutiveMatrix(i, j) = 4 * C10 / DeterminantC * (InvCdotInvC + InvCtimesInvC) +
                                        4 * C01 * (ItimesI - DC_DC - (ItimesInvC + InvCtimesI) / DeterminantC + FirstInvariant / DeterminantC * (InvCdotInvC + InvCtimesInvC)); 
           
        }
    }

   /*  BoundedMatrix<double, 3, 3> eigen_vector_matrix, eigen_values_matrix;
    MathUtils<double>::EigenSystem<3>(rConstitutiveMatrix, eigen_vector_matrix, eigen_values_matrix, 1.0e-16, 200);
    KRATOS_WATCH(eigen_vector_matrix);
    KRATOS_WATCH(eigen_values_matrix); */
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateConstitutiveMatrixKirchhoff(
    Matrix &rConstitutiveMatrix,
    const double DeterminantF,
    const double C10,
    const double C01)
{
    rConstitutiveMatrix.clear();

    /* for (IndexType i = 0; i < 3; ++i)
    {
        const IndexType i0 = this->msIndexVoigt2D3C[i][0];
        const IndexType i1 = this->msIndexVoigt2D3C[i][1];

        for (IndexType j = 0; j < 3; ++j)
        {
            const IndexType j0 = this->msIndexVoigt2D3C[j][0];
            const IndexType j1 = this->msIndexVoigt2D3C[j][1];

            rConstitutiveMatrix(i, j) = 0.0;
        }
    } */
    KRATOS_ERROR << "CalculateConstitutiveMatrixKirchhoff is not implemeted yet." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculatePK2Stress(
    const Matrix &rCTensor,
    const Matrix &rInvCTensor,
    Vector &rStressVector,
    const double DeterminantC,
    const double FirstInvariant,
    const double C10,
    const double C01)
{
    const SizeType dimension = WorkingSpaceDimension();
    Matrix stress_matrix;
    Matrix I = IdentityMatrix(dimension, dimension);

    stress_matrix = 2 * C10 * (I - rInvCTensor / DeterminantC) + 2 * C01 * (I / DeterminantC + FirstInvariant * I - rCTensor - FirstInvariant / DeterminantC * rInvCTensor);
    rStressVector = MathUtils<double>::StressTensorToVector(stress_matrix, GetStrainSize());
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateKirchhoffStress(
    const Matrix &rBTensor,
    Vector &rStressVector,
    const double DeterminantF,
    const double LameLambda,
    const double LameMu)
{
    /* Matrix stress_matrix;

    const SizeType dimension = WorkingSpaceDimension();

    stress_matrix = ZeroMatrix(dimension);

    rStressVector = MathUtils<double>::StressTensorToVector(stress_matrix, rStressVector.size()); */
    KRATOS_ERROR << "CalculateKirchhoffStress is not implemeted yet." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters &rValues,
    Vector &rStrainVector)
{
    // 1.-Compute total deformation gradient
    const Matrix &F = rValues.GetDeformationGradientF();

    // 2.-Compute e = 0.5*(inv(C) - I)
    const Matrix C_tensor = prod(trans(F), F);
    ConstitutiveLawUtilities<VoigtSize>::CalculateGreenLagrangianStrain(C_tensor, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D::CalculateAlmansiStrain(
    ConstitutiveLaw::Parameters &rValues,
    Vector &rStrainVector)
{
    // 1.-Compute total deformation gradient
    const Matrix &F = rValues.GetDeformationGradientF();

    // 2.-COmpute e = 0.5*(1-inv(B))
    const Matrix B_tensor = prod(F, trans(F));
    ConstitutiveLawUtilities<VoigtSize>::CalculateAlmansiStrain(B_tensor, rStrainVector);
}

} // Namespace Kratos
