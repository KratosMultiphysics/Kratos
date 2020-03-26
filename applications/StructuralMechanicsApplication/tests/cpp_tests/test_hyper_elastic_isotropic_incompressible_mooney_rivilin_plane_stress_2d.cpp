// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

// System includes

// External includes

// Project includes
#include <math.h>
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"
// Application includes

// Constitutive law
#include "custom_constitutive/hyper_elastic_isotropic_incompressible_mooney_rivilin_plane_stress_2d.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

/**
    * Check the correct calculation of the integrated stress with the CL's
    */
KRATOS_TEST_CASE_IN_SUITE(HyperElasticMooneyRivlin, KratosStructuralMechanicsFastSuite)
{

    Model current_model;
    ModelPart &test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Triangle3D3<NodeType> Geom = Triangle3D3<NodeType>(p_node_1, p_node_2, p_node_3);

    Matrix deformation_gradient_ref = IdentityMatrix(3, 3);
    deformation_gradient_ref(0, 0) = 1.5;
    deformation_gradient_ref(0, 1) = 0.5;
    deformation_gradient_ref(1, 0) = 0.2;
    deformation_gradient_ref(1, 1) = 1.2;

    Properties hyp_material_properties;
    ProcessInfo hyp_process_info;
    ConstitutiveLaw::Parameters hyp_cl_parameters;
    Vector hyp_stress_vector, hyp_strain_vector;
    hyp_stress_vector = ZeroVector(3);
    hyp_strain_vector = ZeroVector(3);
    Matrix hyp_deformation_gradient = deformation_gradient_ref;
    Matrix hyp_const_matrix = IdentityMatrix(3, 3);
    Matrix hyp_cauchy_green = prod(trans(hyp_deformation_gradient), hyp_deformation_gradient);
    double C1 = 1.925E5;
    double C2 = 1.925E4;

    double C11 = hyp_cauchy_green(0, 0);
    double C12 = hyp_cauchy_green(0, 1);
    double C22 = hyp_cauchy_green(1, 1);
    double C33 = 1 / (C11 * C22 - C12 * C12);

    hyp_material_properties.SetValue(YOUNG_MODULUS, C1);
    hyp_material_properties.SetValue(POISSON_RATIO, C2);
    auto p_hyp_cl = KratosComponents<ConstitutiveLaw>().Get("HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D").Clone();
    Flags hyp_cl_options;
    hyp_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    hyp_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    hyp_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    hyp_cl_parameters.SetElementGeometry(Geom);
    hyp_cl_parameters.SetMaterialProperties(hyp_material_properties);
    hyp_cl_parameters.SetStressVector(hyp_stress_vector);
    hyp_cl_parameters.SetStrainVector(hyp_strain_vector);
    hyp_cl_parameters.SetDeformationGradientF(hyp_deformation_gradient);
    hyp_cl_parameters.SetProcessInfo(hyp_process_info);
    hyp_cl_parameters.SetOptions(hyp_cl_options);
    hyp_cl_parameters.SetConstitutiveMatrix(hyp_const_matrix);

    Vector hyp_plane_stress_ref = ZeroVector(3);
    //https://ocw.mit.edu/resources/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/nonlinear/lecture-15/MITRES2_002S10_lec15.pdf (Slide 31)
    hyp_plane_stress_ref[0] = 2 * C1 * (1 - C33 * C33 * C22) + 2 * C2 * (C33 * 1 + (1 - C33 * C33 * (C11 + C22)) * C22);
    hyp_plane_stress_ref[1] = 2 * C1 * (1 - C33 * C33 * C11) + 2 * C2 * (C33 * 1 + (1 - C33 * C33 * (C11 + C22)) * C11);
    hyp_plane_stress_ref[2] = 2 * C1 * (0 + C33 * C33 * C12) + 2 * C2 * (C33 * 0 - (1 - C33 * C33 * (C11 + C22)) * C12);

    p_hyp_cl->CalculateMaterialResponsePK2(hyp_cl_parameters);

    Matrix hyp_const_matrix_ref = ZeroMatrix(3, 3);
    // Cijkl = 2*DSij/DCkl from mathematica
    hyp_const_matrix_ref(0, 0) = -8 * C22 * (pow(C12, 2) * C2 + C1 * C22 + C2 * pow(C22, 2)) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(0, 1) = 4 * (-C1 * pow(C12, 2) - 2 * C11 * pow(C12, 2) * C2 + pow(C12, 6) * C2 - C1 * C11 * C22 - 2 * pow(C12, 2) * C2 * C22 - 3 * C11 * pow(C12, 4) * C2 * C22 + 3 * pow(C11, 2) * pow(C12, 2) * C2 * pow(C22, 2) - pow(C11, 3) * C2 * pow(C22, 3)) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(0, 2) = 4 * C12 * (pow(C12, 2) * C2 + 2 * C1 * C22 + C11 * C2 * C22 + 2 * C2 * pow(C22, 2)) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(1, 1) = -8 * C11 * (C1 * C11 + pow(C11, 2) * C2 + pow(C12, 2) * C2) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(1, 2) = 4 * C12 * (2 * C1 * C11 + 2 * pow(C11, 2) * C2 + pow(C12, 2) * C2 + C11 * C2 * C22) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(2, 2) = -2 * (3 * C1 * pow(C12, 2) + 3 * C11 * pow(C12, 2) * C2 + pow(C12, 6) * C2 + C1 * C11 * C22 + pow(C11, 2) * C2 * C22 + 3 * pow(C12, 2) * C2 * C22 - 3 * C11 * pow(C12, 4) * C2 * C22 + C11 * C2 * pow(C22, 2) + 3 * pow(C11, 2) * pow(C12, 2) * C2 * pow(C22, 2) - pow(C11, 3) * C2 * pow(C22, 3)) / pow((pow(C12, 2) - C11 * C22), 3);
    hyp_const_matrix_ref(1, 0) = hyp_const_matrix_ref(0, 1);
    hyp_const_matrix_ref(2, 0) = hyp_const_matrix_ref(0, 2);
    hyp_const_matrix_ref(2, 1) = hyp_const_matrix_ref(1, 2);

    std::cout << "################ Final Results ##################" << std::endl;
    KRATOS_WATCH(hyp_strain_vector);
    KRATOS_WATCH(hyp_stress_vector);
    KRATOS_WATCH(hyp_plane_stress_ref);
    KRATOS_WATCH(hyp_const_matrix_ref);
    KRATOS_WATCH(hyp_const_matrix);

    /* double stretch = hyp_deformation_gradient(0, 0);
    double equi_biaxial_stress = 2 * (pow(stretch, 2) - pow(stretch, -4)) * (C1 + pow(stretch, 2) * C2);
    double equi_biaxial_stress_pk_2 = equi_biaxial_stress / pow(stretch, 2);
    KRATOS_WATCH(equi_biaxial_stress_pk_2); */

    BoundedMatrix<double, 2, 2> hyp_eigen_vector_matrix, hyp_eigen_values_matrix;
    MathUtils<double>::EigenSystem<2>(hyp_cauchy_green, hyp_eigen_vector_matrix, hyp_eigen_values_matrix, 1.0e-16, 200);
    double lambda_1 = std::sqrt(hyp_eigen_values_matrix(0, 0));
    double lambda_2 = std::sqrt(hyp_eigen_values_matrix(1, 1));
    // From  Finite Element Analysis of Solids and Structures - M . A. Crisfield
    double S1_ref = 2 * C1 * (1 - pow(lambda_1, -4) * pow(lambda_2, -2)) + 2 * C2 * (-pow(lambda_1, -4) + pow(lambda_2, 2));
    double S2_ref = 2 * C1 * (1 - pow(lambda_2, -4) * pow(lambda_1, -2)) + 2 * C2 * (-pow(lambda_2, -4) + pow(lambda_1, 2));

    Matrix hyp_stress_matrix = ZeroMatrix(2, 2);
    hyp_stress_matrix(0, 0) = hyp_stress_vector[0];
    hyp_stress_matrix(1, 1) = hyp_stress_vector[1];
    hyp_stress_matrix(0, 1) = hyp_stress_vector[2];
    hyp_stress_matrix(1, 0) = hyp_stress_matrix(0, 1);
    BoundedMatrix<double, 2, 2> stress_eigen_vector_matrix, stress_eigen_values_matrix;
    MathUtils<double>::EigenSystem<2>(hyp_stress_matrix, stress_eigen_vector_matrix, stress_eigen_values_matrix, 1.0e-16, 200);
    double S1 = stress_eigen_values_matrix(0, 0);
    double S2 = stress_eigen_values_matrix(1, 1);
    KRATOS_WATCH(S1_ref);
    KRATOS_WATCH(S1);
    KRATOS_WATCH(S2_ref);
    KRATOS_WATCH(S2);

    double V11 = stress_eigen_vector_matrix(0, 0);
    double V22 = stress_eigen_vector_matrix(1, 1);
    double V12 = stress_eigen_vector_matrix(0, 1);
    double V21 = stress_eigen_vector_matrix(1, 0);

    double D1111 = hyp_const_matrix(0, 0);
    double D1122 = hyp_const_matrix(0, 1);
    double D1112 = hyp_const_matrix(0, 2);
    double D2211 = hyp_const_matrix(1, 0);
    double D2222 = hyp_const_matrix(1, 1);
    double D2212 = hyp_const_matrix(1, 2);
    double D1211 = hyp_const_matrix(2, 0);
    double D1222 = hyp_const_matrix(2, 1);
    double D1212 = hyp_const_matrix(2, 2);

    double D11 = pow(V11, 2) * D1111 * pow(V11, 2) + pow(V11, 2) * D1122 * pow(V12, 2) + 2 * pow(V11, 2) * D1112 * V11 * V12 + pow(V12, 2) * D2211 * pow(V11, 2) + pow(V12, 2) * D2222 * pow(V12, 2) +
                 2 * pow(V12, 2) * D2212 * V11 * V12 + 2 * V11 * V12 * D1211 * pow(V11, 2) + 2 * V11 * V12 * D1222 * pow(V12, 2) + 4 * V11 * V12 * D1212 * V11 * V12;
    double D22 = pow(V21, 2) * D1111 * pow(V21, 2) + pow(V21, 2) * D1122 * pow(V22, 2) + 2 * pow(V21, 2) * D1112 * V21 * V22 + pow(V22, 2) * D2211 * pow(V21, 2) + pow(V22, 2) * D2222 * pow(V22, 2) +
                 2 * pow(V22, 2) * D2212 * V21 * V22 + 2 * V21 * V22 * D1211 * pow(V21, 2) + 2 * V21 * V22 * D1222 * pow(V22, 2) + 4 * V21 * V22 * D1212 * V21 * V22;
    double D12 = pow(V11, 2) * D1111 * pow(V21, 2) + pow(V11, 2) * D1122 * pow(V22, 2) + 2 * pow(V11, 2) * D1112 * V21 * V22 + pow(V12, 2) * D2211 * pow(V21, 2) + pow(V12, 2) * D2222 * pow(V22, 2) +
                 2 * pow(V12, 2) * D2212 * V21 * V22 + 2 * V11 * V12 * D1211 * pow(V21, 2) + 2 * V11 * V12 * D1222 * pow(V22, 2) + 4 * V11 * V12 * D1212 * V21 * V22;
    // From  Finite Element Analysis of Solids and Structures - M . A. Crisfield
    double D11_ref = 8 * pow(lambda_1, -6) * (C1 * pow(lambda_2, -2) + C2);
    double D22_ref = 8 * pow(lambda_2, -6) * (C1 * pow(lambda_1, -2) + C2);
    double D12_ref = 4 * (C1 * pow(lambda_1, -4) * pow(lambda_2, -4) + C2);
    KRATOS_WATCH(stress_eigen_vector_matrix);
    KRATOS_WATCH(hyp_eigen_vector_matrix);
    KRATOS_WATCH(D11);
    KRATOS_WATCH(D11_ref);

    KRATOS_WATCH(D22);
    KRATOS_WATCH(D22_ref);
    KRATOS_WATCH(D12);
    KRATOS_WATCH(D12_ref);
    std::cout << "##################  ################" << std::endl;

    // Perturbation analysis to verify tangent matrix D
    double perturbation_threshold = 1e-8;
    Properties material_properties;
    ProcessInfo process_info;
    ConstitutiveLaw::Parameters cl_parameters;
    Vector stress_vector_unpert, strain_vector_unpert, stress_vector_pert, strain_vector_pert, strain_perturbation;
    Matrix deformation_gradient_unpert = IdentityMatrix(3, 3);
    Matrix deformation_gradient_pert = IdentityMatrix(3, 3);
    Matrix const_matrix = IdentityMatrix(3, 3);
    Matrix cauchy_green_pert, green_lagrange_unpert;
    cauchy_green_pert = IdentityMatrix(3, 3);
    Matrix I = IdentityMatrix(3, 3);
    deformation_gradient_unpert(0, 0) = 1.5;
    deformation_gradient_unpert(0, 1) = 0.5;
    deformation_gradient_unpert(1, 0) = deformation_gradient_unpert(0, 1);
    deformation_gradient_unpert(1, 1) = 1.2;

    material_properties.SetValue(YOUNG_MODULUS, C1);
    material_properties.SetValue(POISSON_RATIO, C2);
    auto p_cl = KratosComponents<ConstitutiveLaw>().Get("HyperElasticIsotropicIncompressibleMooneyRivlinPlaneStress2D").Clone();
    //auto p_cl = KratosComponents<ConstitutiveLaw>().Get("LinearElasticPlaneStress2D").Clone();
    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStressVector(stress_vector_pert);
    cl_parameters.SetStrainVector(strain_vector_pert);
    cl_parameters.SetDeformationGradientF(deformation_gradient_unpert);
    cl_parameters.SetProcessInfo(process_info);
    cl_parameters.SetOptions(cl_options);
    cl_parameters.SetConstitutiveMatrix(const_matrix);
    p_cl->CalculateMaterialResponsePK2(cl_parameters);

    stress_vector_unpert = stress_vector_pert;
    strain_vector_unpert = strain_vector_pert;
    green_lagrange_unpert = 0.5 * (prod(trans(deformation_gradient_unpert), deformation_gradient_unpert) - I);

    KRATOS_WATCH(stress_vector_unpert);
    KRATOS_WATCH(strain_vector_unpert);
    KRATOS_WATCH(green_lagrange_unpert);
    int exp = 1;
    double perturbation = 1.0;
    std::vector<double> error_norm_vector;
    std::vector<double> perturbation_vector;

    // D matrix is calculated only once
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    while (perturbation > perturbation_threshold)
    {
        perturbation = pow(10, -exp);
        strain_perturbation = perturbation * strain_vector_unpert;

        strain_vector_pert = strain_vector_unpert + strain_perturbation;

        cauchy_green_pert(0, 0) = 2 * strain_vector_pert[0] + 1;
        cauchy_green_pert(1, 1) = 2 * strain_vector_pert[1] + 1;
        cauchy_green_pert(0, 1) = strain_vector_pert[2]; // gammma = 2* E
        cauchy_green_pert(1, 0) = cauchy_green_pert(0, 1);
        BoundedMatrix<double, 3, 3> eigen_vector_matrix, eigen_values_matrix;
        MathUtils<double>::EigenSystem<3>(cauchy_green_pert, eigen_vector_matrix, eigen_values_matrix, 1.0e-16, 200);
        for (unsigned int i = 0; i < 3; ++i)
            eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));
        deformation_gradient_pert = prod(trans(eigen_vector_matrix), eigen_values_matrix);
        deformation_gradient_pert = prod(deformation_gradient_pert, eigen_vector_matrix);
        cl_parameters.SetDeformationGradientF(deformation_gradient_pert);
        p_cl->CalculateMaterialResponsePK2(cl_parameters);
        Vector calc_stress_diff = stress_vector_pert - stress_vector_unpert;
        Vector pred_stress_diff = prod(const_matrix, strain_perturbation);
        Vector error_vector = calc_stress_diff - pred_stress_diff;
        double error_norm = norm_2(error_vector);
        perturbation_vector.push_back(perturbation);
        error_norm_vector.push_back(error_norm);
        KRATOS_WATCH(strain_vector_pert);
        KRATOS_WATCH(pred_stress_diff);
        KRATOS_WATCH(calc_stress_diff);
        exp += 1;
    }

    std::cout << "Perturbation\t"
              << "Error" << std::endl;

    for (unsigned int i = 0; i < perturbation_vector.size(); ++i)
        std::cout << perturbation_vector[i] << "\t\t" << error_norm_vector[i] << std::endl;

    // Check the stress results
    for (unsigned int i = 0; i < 3; i++)
    {
        KRATOS_CHECK_NEAR(hyp_stress_vector[i], hyp_plane_stress_ref[i], 0.0001e6);
    }

    // Check the constitutive matrix results
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
            KRATOS_CHECK_NEAR(hyp_const_matrix(i, j), hyp_const_matrix_ref(i, j), 0.0001e6);
    }
}

} // namespace Testing
} // namespace Kratos
