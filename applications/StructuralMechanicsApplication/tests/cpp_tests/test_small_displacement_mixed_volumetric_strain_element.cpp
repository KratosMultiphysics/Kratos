// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
// #include "includes/gid_io.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"
#include "factories/linear_solver_factory.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

namespace Kratos
{

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

namespace Testing
{
    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElement2D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);

        p_element->Initialize(r_model_part.GetProcessInfo()); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({51153.8, 51153.8, -1822.18, -12692.3, -38461.5, 10548.1, -38461.5, -12692.3, 3966.35});
        const std::vector<double> expected_LHS_row_0({778846, 9615.38, -317308, -394231, -384615, -317308, -384615, 375000, -317308});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElement2D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D4N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_4->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;
        p_node_4->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.03;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(12);
        Matrix LHS = ZeroMatrix(12,12);

        p_element->Initialize(r_model_part.GetProcessInfo()); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({-23221.2, 56875, -51086.2, 23221.2, 55288.5, -54338.5, 18461.5, -55288.5, -33158.8, -18461.5, -56875, -18483.7});
        const std::vector<double> expected_LHS_row_0({519231,4807.69,-317308,-134615,-379808,-317308,-259615,-4807.69,-158654,-125000,379808,-158654});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElement3D8N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
        auto p_node_5 = r_model_part.CreateNewNode(5, 0.0 , 0.0 , 1.0);
        auto p_node_6 = r_model_part.CreateNewNode(6, 1.0 , 0.0 , 1.0);
        auto p_node_7 = r_model_part.CreateNewNode(7, 1.0 , 1.0 , 1.0);
        auto p_node_8 = r_model_part.CreateNewNode(8, 0.0 , 1.0 , 1.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4,5,6,7,8};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D8N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_4->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;
        p_node_4->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.03;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(32);
        Matrix LHS = ZeroMatrix(32,32);

        p_element->Initialize(r_model_part.GetProcessInfo()); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({-4144.23, 16003.2, -12609, -7830.25, 4144.23, 15544.9, -12838.1, -9873.46, 2769.23, -34775.6, 7309.29, 2080.25, -2769.23, -35234, 7767.63, 9947.53, -2072.12, 17617, -6621.79, -17148.1, 2072.12, 17387.8, -6392.63, -17224.5, 1384.62, 1842.95, 11921.5, -15028.5, -1384.62, 1613.78, 11463.1, -12985.3});
        const std::vector<double> expected_LHS_row_0({286752,22756.4,22756.4,-91666.7,-30341.9,-105449,-105449,-91666.7,-79273.5,-22756.4,-52724.4,-45833.3,15170.9,105449,11378.2,-45833.3,15170.9,11378.2,105449,-45833.3,-79273.5,-52724.4,-22756.4,-45833.3,-71688,-11378.2,-11378.2,-22916.7,-56517.1,52724.4,52724.4,-22916.7});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    void LinearPerturbationField(
        ModelPart &rModelPart,
        const double c_1,
        const double c_2,
        const double c_3 = 0.0)
    {
        array_1d<double, 3> aux_disp;
        for (auto &r_node : rModelPart.Nodes()) {
            aux_disp[0] = c_1 * r_node.X0();
            aux_disp[1] = c_2 * r_node.Y0();
            aux_disp[2] = c_3 * r_node.Z0();
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = c_1 + c_2 + c_3;

            r_node.X() = (1 + c_1) * r_node.X0();
            r_node.Y() = (1 + c_2) * r_node.Y0();
            r_node.Z() = (1 + c_3) * r_node.Z0();
        }
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * LHS and RHS triangle element test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElement2D3NResidual, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.1 , 0.1 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.5 , 0.2 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.9 , 0.4 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 1, element_nodes, p_elem_prop);

        // Initialize the element to initialize the constitutive law
        p_element->Initialize(r_model_part.GetProcessInfo());

        // Set a fake displacement and volumetric strain field to compute the residual
        const double alpha = -2.0e-5;
        const double beta = 1.0e-5;
        LinearPerturbationField(r_model_part, alpha, beta);

        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Perturb the previous displacement and volumetric strain field to compute the residual
        const double alpha_perturbed = 1.25e-5;
        const double beta_perturbed = 0.25e-5;
        LinearPerturbationField(r_model_part, alpha_perturbed, beta_perturbed);

        Vector RHS_perturbed = ZeroVector(9);
        p_element->CalculateRightHandSide(RHS_perturbed, r_model_part.GetProcessInfo());

        // Calculate the perturbation RHS
        const double delta_alpha = alpha_perturbed - alpha;
        const double delta_beta = beta_perturbed - beta;
        LinearPerturbationField(r_model_part, delta_alpha, delta_beta);

        Vector RHS_delta = ZeroVector(9);
        p_element->CalculateRightHandSide(RHS_delta, r_model_part.GetProcessInfo());

        // Check the error
        const Vector RHS_error = RHS_perturbed - (RHS + RHS_delta);

        // Check the LHS
        array_1d<double, 9> perturbation_vector = ZeroVector(9);
        for (auto &r_node : r_model_part.Nodes()) {
            perturbation_vector[(r_node.Id() - 1) * 3] = delta_alpha * r_node.X0();
            perturbation_vector[(r_node.Id() - 1) * 3 + 1] = delta_beta * r_node.Y0();
            perturbation_vector[(r_node.Id() - 1) * 3 + 2] = delta_alpha + delta_beta;
        }
        const Vector RHS_from_LHS = RHS - prod(LHS, perturbation_vector);
        const Vector RHS_from_LHS_error = RHS_perturbed - RHS_from_LHS;

        // Check RHS and LHS results
        const double tolerance = 1.0e-12;
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_error, ZeroVector(r_model_part.NumberOfNodes() * 3), tolerance);
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_from_LHS_error, ZeroVector(r_model_part.NumberOfNodes() * 3), tolerance);
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * LHS and RHS tetrahedra element test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElement3D4NResidual, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.1 , 0.1 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.5 , 0.2 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.9 , 0.4 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.3 , 0.3 , 0.5);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 1, element_nodes, p_elem_prop);

        // Initialize the element to initialize the constitutive law
        p_element->Initialize(r_model_part.GetProcessInfo());

        // Set a fake displacement and volumetric strain field to compute the residual
        const double alpha = -2.0;
        const double beta = 1.0;
        const double gamma = 0.5;
        LinearPerturbationField(r_model_part, alpha, beta, gamma);

        Vector RHS = ZeroVector(16);
        Matrix LHS = ZeroMatrix(16,16);
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Perturb the previous displacement and volumetric strain field to compute the residual
        const double alpha_perturbed = 1.25;
        const double beta_perturbed = 0.25;
        const double gamma_perturbed = -0.5;
        LinearPerturbationField(r_model_part, alpha_perturbed, beta_perturbed, gamma_perturbed);

        Vector RHS_perturbed = ZeroVector(16);
        p_element->CalculateRightHandSide(RHS_perturbed, r_model_part.GetProcessInfo());

        // Calculate the perturbation RHS
        const double delta_alpha = alpha_perturbed - alpha;
        const double delta_beta = beta_perturbed - beta;
        const double delta_gamma = gamma_perturbed - gamma;
        LinearPerturbationField(r_model_part, delta_alpha, delta_beta, delta_gamma);

        Vector RHS_delta = ZeroVector(16);
        p_element->CalculateRightHandSide(RHS_delta, r_model_part.GetProcessInfo());

        // Check the error
        const Vector RHS_error = RHS_perturbed - (RHS + RHS_delta);

        // Check the LHS
        array_1d<double, 16> perturbation_vector = ZeroVector(16);
        for (auto &r_node : r_model_part.Nodes()) {
            perturbation_vector[(r_node.Id() - 1) * 4] = delta_alpha * r_node.X0();
            perturbation_vector[(r_node.Id() - 1) * 4 + 1] = delta_beta * r_node.Y0();
            perturbation_vector[(r_node.Id() - 1) * 4 + 2] = delta_gamma * r_node.Z0();
            perturbation_vector[(r_node.Id() - 1) * 4 + 3] = delta_alpha + delta_beta + delta_gamma;
        }
        const Vector RHS_from_LHS = RHS - prod(LHS, perturbation_vector);
        const Vector RHS_from_LHS_error = RHS_perturbed - RHS_from_LHS;

        // Check RHS and LHS results
        const double tolerance = 1.0e-9;
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_error, ZeroVector(r_model_part.NumberOfNodes() * 4), tolerance);
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_from_LHS_error, ZeroVector(r_model_part.NumberOfNodes() * 4), tolerance);
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * Zienkiewicz patch test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainElementZienkiewiczPatch, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(REACTION);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(REACTION_STRAIN);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(1);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+02);
        p_elem_prop->SetValue(POISSON_RATIO, 0.4995);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.0 , 1.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 1.0 , 0.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes_1 {1,3,2};
        std::vector<ModelPart::IndexType> element_nodes_2 {1,4,3};
        auto p_element_1 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 1, element_nodes_1, p_elem_prop);
        auto p_element_2 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 2, element_nodes_2, p_elem_prop);

        // Create the load condition
        auto p_cond_prop = r_model_part.CreateNewProperties(0);
        std::vector<ModelPart::IndexType> condition_nodes_1 {3};
        auto p_condition_1 = r_model_part.CreateNewCondition("PointLoadCondition2D1N", 1, condition_nodes_1, p_cond_prop);
        array_1d<double,3> point_load = ZeroVector(3);
        point_load[1] = 1.0;
        p_condition_1->SetValue(POINT_LOAD, point_load);

        // Add DOFs
        for (auto &r_node: r_model_part.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            r_node.AddDof(VOLUMETRIC_STRAIN, REACTION_STRAIN);
        }

        // Initialize the elements to initialize the constitutive law
        for (auto &r_elem : r_model_part.Elements()) {
            r_elem.Initialize(r_model_part.GetProcessInfo());
        }

        // Construct the linear solver pointer
        Parameters linear_solver_settings(R"({"solver_type": "skyline_lu_factorization"})");
        LinearSolverFactory<SparseSpaceType, LocalSpaceType> linear_solver_factory;
        auto p_linear_solver = linear_solver_factory.Create(linear_solver_settings);

        // Create the linear strategy
        Scheme<SparseSpaceType,LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>>();

        bool calculate_norm_dx = false;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(p_linear_solver);

        auto p_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            r_model_part,
            p_scheme,
            p_linear_solver,
            p_builder_and_solver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx);

        p_solving_strategy->Check();

        // Fix the boundary nodes
        p_node_1->Fix(DISPLACEMENT_X);
        p_node_1->Fix(DISPLACEMENT_Y);
        p_node_2->Fix(DISPLACEMENT_X);
        p_node_2->Fix(DISPLACEMENT_Y);
        p_node_4->Fix(DISPLACEMENT_X);
        p_node_4->Fix(DISPLACEMENT_Y);

        // Solve the problem
        p_solving_strategy->Solve();

        // Check results
        const double tolerance = 1.0e-6;
        const double expected_vol_strain = 1.4965e-05;
        const std::vector<double> expected_disp = {-0.000318007, 0.000347937, 0};
        KRATOS_CHECK_VECTOR_NEAR(p_node_3->FastGetSolutionStepValue(DISPLACEMENT), expected_disp, tolerance);
        KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN), expected_vol_strain, tolerance);

        // // GiD output
        // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/test_small_displacement_mixed_volumetric_strain_element_zienkiewicz_patch", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		// gid_io_fluid.InitializeMesh(0.0);
		// gid_io_fluid.WriteMesh(r_model_part.GetMesh());
		// gid_io_fluid.FinalizeMesh();
		// gid_io_fluid.InitializeResults(0, r_model_part.GetMesh());
		// gid_io_fluid.WriteNodalResults(DISPLACEMENT, r_model_part.Nodes(), 0, 0);
		// gid_io_fluid.WriteNodalResults(VOLUMETRIC_STRAIN, r_model_part.Nodes(), 0, 0);
		// gid_io_fluid.FinalizeResults();

    }

} // namespace Testing
} // namespace Kratos.
