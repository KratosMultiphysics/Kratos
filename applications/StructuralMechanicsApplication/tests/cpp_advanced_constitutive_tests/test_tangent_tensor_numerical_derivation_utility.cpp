// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{
namespace Testing
{

/// Node definition
typedef Node<3> NodeType;

/**
 * @brief This sets the basic case
 */
void SettingBasicCase(
    ModelPart& rModelPart,
    ConstitutiveLaw::Parameters& rCLConfigurationValues,
    Properties& rMaterialProperties,
    Vector& rStressVector,
    Vector& rStrainVector,
    Matrix& rTangentModuli,
    Matrix& rF,
    double& rDetF,
    const bool ElementProvidedStrain = true,
    const std::size_t Case = 1
    )
{
    ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    Flags constitutive_law_options;
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, ElementProvidedStrain);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    rMaterialProperties.SetValue(YOUNG_MODULUS, 210e9);
    rMaterialProperties.SetValue(POISSON_RATIO, 0.22);

    rStressVector = ZeroVector(6);

    rStrainVector = ZeroVector(6); // Case 0
    if (Case == 1) {
        rStrainVector[2] = 8.0e-5;
        rStrainVector[5] = 1.6941e-21;
    } else if (Case == 2) {
        rMaterialProperties.SetValue(YOUNG_MODULUS, 10e9);
        rMaterialProperties.SetValue(POISSON_RATIO, 0.0);
        rMaterialProperties.SetValue(YIELD_STRESS, 1.0e5);
        rMaterialProperties.SetValue(FRICTION_ANGLE, 32.0);
        rMaterialProperties.SetValue(DILATANCY_ANGLE, 16.0);
        rMaterialProperties.SetValue(SOFTENING_TYPE, 1);
        rMaterialProperties.SetValue(FRACTURE_ENERGY, 100.0);
        rMaterialProperties.SetValue(HARDENING_CURVE, 3);

        rStrainVector[1] = 1.0e-5;
    }

    // Compute equivalent F
    rF = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i) {
        rF(i, i) = 1.0 + rStrainVector[i];
    }

    for (int i = 3; i < 6; ++i) {
        const int equivalent_i = (i == 3) ? 0 : (i == 4) ? 1 : 0;
        const int equivalent_j = (i == 3) ? 1 : 2;
        rF(equivalent_i, equivalent_j) = 0.5 * rStrainVector[i];
        rF(equivalent_j, equivalent_i) = 0.5 * rStrainVector[i];
    }
    rDetF = MathUtils<double>::DetMat(rF);

    rTangentModuli = ZeroMatrix(6, 6);

    rCLConfigurationValues.SetMaterialProperties(rMaterialProperties);
    rCLConfigurationValues.SetDeformationGradientF(rF);
    rCLConfigurationValues.SetDeterminantF(rDetF);
    rCLConfigurationValues.SetStrainVector(rStrainVector);
    rCLConfigurationValues.SetStressVector(rStressVector);
    rCLConfigurationValues.SetOptions(constitutive_law_options);
    rCLConfigurationValues.SetProcessInfo(r_current_process_info);
    rCLConfigurationValues.SetConstitutiveMatrix(rTangentModuli);
}

/**
 * @brief This method creates a simple geometry (hexahedra)
 */
void Create3DGeometryHexahedra(
    ModelPart& rThisModelPart,
    const std::string ElementName = "SmallDisplacementElement3D8N"
    )
{
    auto& process_info = rThisModelPart.GetProcessInfo();
    process_info[STEP] = 2;
    process_info[NL_ITERATION_NUMBER] = 2;

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);

    // Now we create the elements
    rThisModelPart.CreateNewElement(ElementName, 1, {{5,8,6,2,3,7,4,1}}, p_elem_prop);

    // Initialize elements
    for (auto& elem : rThisModelPart.Elements()) {
        elem.Initialize();
        elem.InitializeSolutionStep(process_info);
        elem.InitializeNonLinearIteration(process_info);
    }
}

/**
 * @brief This method creates a simple geometry (tetrahedra)
 */
void Create3DGeometryTetrahedra(
    ModelPart& rThisModelPart,
    const std::string ElementName = "SmallDisplacementElement3D4N"
    )
{
    auto& process_info = rThisModelPart.GetProcessInfo();
    process_info[STEP] = 2;
    process_info[NL_ITERATION_NUMBER] = 2;

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);
    // First we create the nodes
    rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    // Now we create the elements
    rThisModelPart.CreateNewElement(ElementName, 1, {{12,10,8,9}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 2, {{4,6,9,7}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 3, {{11,7,9,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 4, {{5,3,8,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 5, {{4,6,7,3}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 6, {{2,3,5,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 7, {{10,9,6,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 8, {{7,8,3,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 9, {{7,8,6,9}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 10, {{4,1,6,3}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 11, {{9,12,11,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 12, {{3,2,1,6}}, p_elem_prop);

    // Initialize elements
    for (auto& elem : rThisModelPart.Elements()) {
        elem.Initialize();
        elem.InitializeSolutionStep(process_info);
        elem.InitializeNonLinearIteration(process_info);
    }
}

/**
 * @brief This computes the convergence rate of the CL
 */
template<std::size_t TDim = 3>
void ComputingConvergenceRate(
    ConstitutiveLaw::Pointer pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rCLConfigurationValues,
    Vector& rStressVector,
    Vector& rStrainVector,
    Matrix& rTangentModuli,
    Matrix& rF,
    double& rDetF,
    const bool FiniteDeformation = false,
    const bool Debug = false,
    const bool AnalyticReference = true,
    double Alpha = 1.0,
    const std::size_t MaxNumberIters = 10,
    const std::size_t IntegrationOrder = 2
    )
{
    const std::size_t voigt_size = TDim == 3 ? 6 : 3;

    // The delta parameters
    const Vector initial_strain_vector = rStrainVector;
    const Vector initial_stress_vector = rStressVector;
    const Matrix initial_deformation_gradient_F = rF;
    Vector delta_strain_vector = rStrainVector;
    Matrix delta_deformation_gradient_F = rF;

    // Ensure the proper flag
    Flags& cl_options = rCLConfigurationValues.GetOptions();
    const bool use_element_provided_strain = cl_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // First error computation
    Vector expected_delta_stress = ZeroVector(voigt_size);
    Vector computed_delta_stress = ZeroVector(voigt_size);

    const double quadratic_threshold = IntegrationOrder == 2 ? 1.8 : 1.0;

    std::vector<double> vector_errors(MaxNumberIters);
    for (std::size_t iter = 0; iter < MaxNumberIters; ++iter) {
        Alpha *= 0.5;
        noalias(delta_strain_vector) = Alpha * initial_strain_vector;
        noalias(delta_deformation_gradient_F) = Alpha * (initial_deformation_gradient_F - IdentityMatrix(TDim)) + IdentityMatrix(TDim);
        noalias(rStrainVector) =  initial_strain_vector + delta_strain_vector;
        if (!use_element_provided_strain) {
            noalias(rF) = prod(delta_deformation_gradient_F, initial_deformation_gradient_F);
            rDetF = MathUtils<double>::DetMat(rF);
            rCLConfigurationValues.SetDeformationGradientF(rF);
            rCLConfigurationValues.SetDeterminantF(rDetF);
        }

        pConstitutiveLaw->CalculateMaterialResponse(rCLConfigurationValues, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
        noalias(expected_delta_stress) = (rStressVector - initial_stress_vector);
        if (!use_element_provided_strain) {
            noalias(delta_strain_vector) = rStrainVector - initial_strain_vector;
        }

        if (FiniteDeformation) {
            TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(rCLConfigurationValues, pConstitutiveLaw.get(), ConstitutiveLaw::StressMeasure::StressMeasure_PK2, true, IntegrationOrder);
        } else {
            TangentOperatorCalculatorUtility::CalculateTangentTensor(rCLConfigurationValues, pConstitutiveLaw.get(), ConstitutiveLaw::StressMeasure::StressMeasure_PK2, true, IntegrationOrder);
        }

        const Matrix tangent_moduli_perturbed = rTangentModuli;

        double error;
        if (AnalyticReference) {
            cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
            pConstitutiveLaw->CalculateMaterialResponse(rCLConfigurationValues, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
            cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            const Matrix aux_error_matrix = tangent_moduli_perturbed - rTangentModuli;
            error = norm_frobenius(aux_error_matrix)/norm_frobenius(rTangentModuli);
        } else {
            noalias(computed_delta_stress) = prod(rTangentModuli, delta_strain_vector);
            const Vector aux_error_vector = computed_delta_stress - expected_delta_stress;
            error = norm_2(aux_error_vector);
        }

        vector_errors[iter] = error;

        if (Debug) {
            KRATOS_WATCH(error)
            if (pConstitutiveLaw->Has(PLASTIC_DISSIPATION)) {
                auto p_clone = pConstitutiveLaw->Clone();
                double plastic_dissipation;
                p_clone->FinalizeMaterialResponse(rCLConfigurationValues, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
                p_clone->GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
                double uniaxial_stress;
                p_clone->CalculateValue(rCLConfigurationValues, UNIAXIAL_STRESS, uniaxial_stress);
                KRATOS_WATCH(plastic_dissipation)
                KRATOS_WARNING_IF("PLASTICITY", plastic_dissipation < 1.0e-12) << "This test is not in plastic range" << std::endl;
            }
        }
    }

    if (MaxNumberIters > 5) {
        for (std::size_t i = 1; i < MaxNumberIters - 4; ++i) { // We discard the first solution
            if ((vector_errors[i + 3] + vector_errors[i + 2] + vector_errors[i + 1] + vector_errors[i]) > 1.0e-4) { // If zero means exact solution
                const double slope = std::log((vector_errors[i + 3] - vector_errors[i + 2])/(vector_errors[i + 2] - vector_errors[i + 1]))/std::log((vector_errors[i + 2] - vector_errors[i + 1])/(vector_errors[i + 1] - vector_errors[i + 0]));

                if (Debug) {
                    KRATOS_WATCH(slope)
                } else { // Check
                    KRATOS_CHECK_GREATER_EQUAL(slope, quadratic_threshold);
                }
            }
        }
    } else {
        const double log_coeff = std::log(0.5);
        for (std::size_t i = 1; i < MaxNumberIters - 1; ++i) { // We discard the first solution
            const double slope = std::log(vector_errors[i]/vector_errors[i-1])/log_coeff;
            if (Debug) {
                KRATOS_WATCH(slope)
            } else { // Check
                KRATOS_CHECK_GREATER_EQUAL(slope, quadratic_threshold);
            }
        }
    }
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the elastic linear tensor
 */
KRATOS_TEST_CASE_IN_SUITE(LinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main");

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_Cauchy);

    Matrix C = ZeroMatrix(6, 6);
    C = p_constitutive_law->CalculateValue(cl_configuration_values,CONSTITUTIVE_MATRIX, C);

    TangentOperatorCalculatorUtility::CalculateTangentTensor(cl_configuration_values, p_constitutive_law.get());

    const double tolerance = 1.0e-4;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (std::abs(C(i, j)) > 0.0) {
                KRATOS_CHECK_LESS_EQUAL((tangent_moduli(i, j) - C(i, j))/C(i, j), tolerance);
            } else if (std::abs(tangent_moduli(i, j)) > 1.0e-6) {
                KRATOS_WARNING("LinearElasticCasePertubationTensorUtility") << "Be careful tangent_moduli(" << i << " ," << j << ") is greater than 0: " <<  tangent_moduli(i, j) << std::endl;
            }
        }
    }
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the hyperelastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(HyperElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main");

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("KirchhoffSaintVenant3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    Matrix C = ZeroMatrix(6, 6);
    C = p_constitutive_law->CalculateValue(cl_configuration_values,CONSTITUTIVE_MATRIX_PK2, C);

    TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(cl_configuration_values, p_constitutive_law.get());

    const double tolerance = 1.0e-4;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (std::abs(C(i, j)) > 0.0) {
                KRATOS_CHECK_LESS_EQUAL((tangent_moduli(i, j) - C(i, j))/C(i, j), tolerance);
            } else if (std::abs(tangent_moduli(i, j)) > 1.0e-6) {
                KRATOS_WARNING("HyperElasticCasePertubationTensorUtility") << "Be careful tangent_moduli(" << i << " ," << j << ") is greater than 0: " <<  tangent_moduli(i, j) << std::endl;
            }
        }
    }
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the elastic linear tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticLinearElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main");

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("LinearElastic3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false, false);
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the Kirchhoff hyperelastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticKirchhoffHyperElasticCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main");

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false);

    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("KirchhoffSaintVenant3DLaw").Clone();
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, true, false);
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the ehyper lastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticSmallStrainIsotropicPlasticity3DVonMisesVonMisesFirstOrderCasePertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, true, 2);

    // Creating constitutive law
    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("SmallStrainIsotropicPlasticity3DVonMisesVonMises").Clone();
    p_material_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);

    // Creating geometry
    Create3DGeometryHexahedra(r_model_part);

    // Assigning geometry
    auto& r_geom = r_model_part.Elements().begin()->GetGeometry();
    cl_configuration_values.SetElementGeometry(r_geom);

    // Initializing law
    p_constitutive_law->InitializeMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    // Computing convergence rate
    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false, false, false, 1.0, 4, 1);
}

/**
 * @brief This test tests that the perturbation utility is valid for computing the ehyper lastic tensor
 */
KRATOS_TEST_CASE_IN_SUITE(QuadraticSmallStrainIsotropicPlasticity3DVonMisesVonMisesCaseSecondOrderPertubationTensorUtility, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    ConstitutiveLaw::Parameters cl_configuration_values;
    Properties::Pointer p_material_properties = r_model_part.CreateNewProperties(1);
    Vector stress_vector, strain_vector;
    Matrix tangent_moduli, deformation_gradient_F;
    double det_deformation_gradient_F;
    SettingBasicCase(r_model_part, cl_configuration_values, *p_material_properties, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, true, 2);

    // Creating constitutive law
    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get("SmallStrainIsotropicPlasticity3DVonMisesVonMises").Clone();
    p_material_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);

    // Creating geometry
    Create3DGeometryHexahedra(r_model_part);

    // Assigning geometry
    auto& r_geom = r_model_part.Elements().begin()->GetGeometry();
    cl_configuration_values.SetElementGeometry(r_geom);

    // Initializing law
    p_constitutive_law->InitializeMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);
    p_constitutive_law->CalculateMaterialResponse(cl_configuration_values, ConstitutiveLaw::StressMeasure::StressMeasure_PK2);

    // Computing convergence rate
    ComputingConvergenceRate(p_constitutive_law, cl_configuration_values, stress_vector, strain_vector, tangent_moduli, deformation_gradient_F, det_deformation_gradient_F, false, false, false, 1.0, 4);
}

} // namespace Testing
} // namespace Kratos
