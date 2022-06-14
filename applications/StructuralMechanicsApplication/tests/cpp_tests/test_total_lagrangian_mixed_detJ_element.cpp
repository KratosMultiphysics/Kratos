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
#include "factories/linear_solver_factory.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/gid_io.h"
#include "processes/structured_mesh_generator_process.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/total_lagrangian_mixed_detJ_element.h"

namespace Kratos
{

namespace Testing
{

// namespace
// {
//     void LinearPerturbationField(
//         ModelPart &rModelPart,
//         const double c_1,
//         const double c_2,
//         const double c_3 = 0.0)
//     {
//         array_1d<double, 3> aux_disp;
//         for (auto &r_node : rModelPart.Nodes()) {
//             aux_disp[0] = c_1 * r_node.X0();
//             aux_disp[1] = c_2 * r_node.Y0();
//             aux_disp[2] = c_3 * r_node.Z0();
//             r_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
//             r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = (1.0+c_1)*(1+c_2)*(1+c_3);
//             r_node.X() = (1 + c_1) * r_node.X0();
//             r_node.Y() = (1 + c_2) * r_node.Y0();
//             r_node.Z() = (1 + c_3) * r_node.Z0();
//         }
//     }
// }

    /**
    * Checks the Total Lagrangian mixed Jacobian determinant element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedDetJElement2D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw"); //FIXME: Use kirchoff in here
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("TotalLagrangianMixedDetJElement2D3N", 1, element_nodes, p_elem_prop);

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

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({58615.4,60579.5,-0.00446042,-19641,-40938.5,0.00879375,-38974.4,-19641,0.00233333});
        const std::vector<double> expected_LHS_row_0({949041,129590,-335577,-402564,-526154,-335577,-546477,396564,-335577});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    // void LinearPerturbationField(
    //     ModelPart &rModelPart,
    //     const double c_1,
    //     const double c_2,
    //     const double c_3 = 0.0)
    // {
    //     array_1d<double, 3> aux_disp;
    //     for (auto &r_node : rModelPart.Nodes()) {
    //         aux_disp[0] = c_1 * r_node.X0();
    //         aux_disp[1] = c_2 * r_node.Y0();
    //         aux_disp[2] = c_3 * r_node.Z0();
    //         r_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    //         r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = c_1 + c_2 + c_3;

    //         r_node.X() = (1 + c_1) * r_node.X0();
    //         r_node.Y() = (1 + c_2) * r_node.Y0();
    //         r_node.Z() = (1 + c_3) * r_node.Z0();
    //     }
    // }

    // /**
    // * Checks the Small Displacement Mixed Strain Element
    // * LHS and RHS triangle element test
    // */
    // KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedDetJElement2D3NResidual, KratosStructuralMechanicsFastSuite)
    // {
    //     Model current_model;
    //     auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

    //     const auto& r_process_info = r_model_part.GetProcessInfo();

    //     r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    //     r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

    //     // Set the element properties
    //     auto p_elem_prop = r_model_part.CreateNewProperties(0);
    //     p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
    //     p_elem_prop->SetValue(POISSON_RATIO, 0.3);
    //     const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw"); //FIXME: use kirchoff here
    //     p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    //     // Create the test element
    //     auto p_node_1 = r_model_part.CreateNewNode(1, 0.1 , 0.1 , 0.0);
    //     auto p_node_2 = r_model_part.CreateNewNode(2, 0.5 , 0.2 , 0.0);
    //     auto p_node_3 = r_model_part.CreateNewNode(3, 0.9 , 0.4 , 0.0);
    //     std::vector<ModelPart::IndexType> element_nodes {1,2,3};
    //     auto p_element = r_model_part.CreateNewElement("TotalLagrangianMixedDetJElement2D3N", 1, element_nodes, p_elem_prop);

    //     // Initialize the element to initialize the constitutive law
    //     p_element->Initialize(r_process_info);

    //     // Set a fake displacement and volumetric strain field to compute the residual
    //     const double alpha = -2.0e-5;
    //     const double beta = 1.0e-5;
    //     LinearPerturbationField(r_model_part, alpha, beta);

    //     Vector RHS = ZeroVector(9);
    //     Matrix LHS = ZeroMatrix(9,9);
    //     p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    //     // Perturb the previous displacement and volumetric strain field to compute the residual
    //     const double alpha_perturbed = 1.25e-5;
    //     const double beta_perturbed = 0.25e-5;
    //     LinearPerturbationField(r_model_part, alpha_perturbed, beta_perturbed);

    //     Vector RHS_perturbed = ZeroVector(9);
    //     p_element->CalculateRightHandSide(RHS_perturbed, r_process_info);

    //     // Calculate the perturbation RHS
    //     const double delta_alpha = alpha_perturbed - alpha;
    //     const double delta_beta = beta_perturbed - beta;
    //     LinearPerturbationField(r_model_part, delta_alpha, delta_beta);

    //     Vector RHS_delta = ZeroVector(9);
    //     p_element->CalculateRightHandSide(RHS_delta, r_process_info);

    //     // Check the error
    //     const Vector RHS_error = RHS_perturbed - (RHS + RHS_delta);

    //     // Check the LHS
    //     array_1d<double, 9> perturbation_vector = ZeroVector(9);
    //     for (auto &r_node : r_model_part.Nodes()) {
    //         perturbation_vector[(r_node.Id() - 1) * 3] = delta_alpha * r_node.X0();
    //         perturbation_vector[(r_node.Id() - 1) * 3 + 1] = delta_beta * r_node.Y0();
    //         perturbation_vector[(r_node.Id() - 1) * 3 + 2] = delta_alpha + delta_beta;
    //     }
    //     const Vector RHS_from_LHS = RHS - prod(LHS, perturbation_vector);
    //     const Vector RHS_from_LHS_error = RHS_perturbed - RHS_from_LHS;

    //     // Check RHS and LHS results
    //     const double tolerance = 1.0e-12;
    //     KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_error, ZeroVector(r_model_part.NumberOfNodes() * 3), tolerance);
    //     KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_from_LHS_error, ZeroVector(r_model_part.NumberOfNodes() * 3), tolerance);
    // }

    // /**
    // * Checks the Small Displacement Mixed Strain Element
    // * LHS and RHS tetrahedra element test
    // */
    // KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedDetJElement3D4NResidual, KratosStructuralMechanicsFastSuite)
    // {
    //     Model current_model;
    //     auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

    //     const auto& r_process_info = r_model_part.GetProcessInfo();

    //     r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    //     r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

    //     // Set the element properties
    //     auto p_elem_prop = r_model_part.CreateNewProperties(0);
    //     p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
    //     p_elem_prop->SetValue(POISSON_RATIO, 0.3);
    //     const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
    //     p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    //     // Create the test element
    //     auto p_node_1 = r_model_part.CreateNewNode(1, 0.1 , 0.1 , 0.0);
    //     auto p_node_2 = r_model_part.CreateNewNode(2, 0.5 , 0.2 , 0.0);
    //     auto p_node_3 = r_model_part.CreateNewNode(3, 0.9 , 0.4 , 0.0);
    //     auto p_node_4 = r_model_part.CreateNewNode(4, 0.3 , 0.3 , 0.5);
    //     std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    //     auto p_element = r_model_part.CreateNewElement("TotalLagrangianMixedDetJElement3D4N", 1, element_nodes, p_elem_prop);

    //     // Initialize the element to initialize the constitutive law
    //     p_element->Initialize(r_process_info);

    //     // Set a fake displacement and volumetric strain field to compute the residual
    //     const double alpha = -2.0;
    //     const double beta = 1.0;
    //     const double gamma = 0.5;
    //     LinearPerturbationField(r_model_part, alpha, beta, gamma);

    //     Vector RHS = ZeroVector(16);
    //     Matrix LHS = ZeroMatrix(16,16);
    //     p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    //     // Perturb the previous displacement and volumetric strain field to compute the residual
    //     const double alpha_perturbed = 1.25;
    //     const double beta_perturbed = 0.25;
    //     const double gamma_perturbed = -0.5;
    //     LinearPerturbationField(r_model_part, alpha_perturbed, beta_perturbed, gamma_perturbed);

    //     Vector RHS_perturbed = ZeroVector(16);
    //     p_element->CalculateRightHandSide(RHS_perturbed, r_process_info);

    //     // Calculate the perturbation RHS
    //     const double delta_alpha = alpha_perturbed - alpha;
    //     const double delta_beta = beta_perturbed - beta;
    //     const double delta_gamma = gamma_perturbed - gamma;
    //     LinearPerturbationField(r_model_part, delta_alpha, delta_beta, delta_gamma);

    //     Vector RHS_delta = ZeroVector(16);
    //     p_element->CalculateRightHandSide(RHS_delta, r_process_info);

    //     // Check the error
    //     const Vector RHS_error = RHS_perturbed - (RHS + RHS_delta);

    //     // Check the LHS
    //     array_1d<double, 16> perturbation_vector = ZeroVector(16);
    //     for (auto &r_node : r_model_part.Nodes()) {
    //         perturbation_vector[(r_node.Id() - 1) * 4] = delta_alpha * r_node.X0();
    //         perturbation_vector[(r_node.Id() - 1) * 4 + 1] = delta_beta * r_node.Y0();
    //         perturbation_vector[(r_node.Id() - 1) * 4 + 2] = delta_gamma * r_node.Z0();
    //         perturbation_vector[(r_node.Id() - 1) * 4 + 3] = delta_alpha + delta_beta + delta_gamma;
    //     }
    //     const Vector RHS_from_LHS = RHS - prod(LHS, perturbation_vector);
    //     const Vector RHS_from_LHS_error = RHS_perturbed - RHS_from_LHS;

    //     // Check RHS and LHS results
    //     const double tolerance = 1.0e-9;
    //     KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_error, ZeroVector(r_model_part.NumberOfNodes() * 4), tolerance);
    //     KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS_from_LHS_error, ZeroVector(r_model_part.NumberOfNodes() * 4), tolerance);
    // }

    /**
     * This test implements the simple patch test example in Bonet and Wood book
     */
    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedDetJElementBonetPatch, KratosStructuralMechanicsFastSuite)
    {
        // Set the test model part
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        r_model_part.AddNodalSolutionStepVariable(REACTION);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(REACTION_STRAIN);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Create the test mesh
        const double min_x = 0.0;
        const double max_x = 1.0;
        const double min_y = 0.0;
        const double max_y = 1.0;
		auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, min_x, min_y, 0.0);
		auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, min_x, max_y, 0.0);
		auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, max_x, max_y, 0.0);
		auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, max_x, min_y, 0.0);
		Quadrilateral2D4<Node<3>> square_geometry(p_point_1, p_point_2, p_point_3, p_point_4);
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions" : 2,
			"element_name" : "TotalLagrangianMixedDetJElement2D3N"
		})");
		StructuredMeshGeneratorProcess(square_geometry, r_model_part, mesher_parameters).Execute();

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(1);
        p_elem_prop->SetValue(YOUNG_MODULUS, 250.0);
        p_elem_prop->SetValue(POISSON_RATIO, 0.25);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("KirchhoffSaintVenantPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());
        for (auto& r_elem: r_model_part.Elements()) {
            r_elem.SetProperties(p_elem_prop);
        }

        // Add DOFs
        for (auto &r_node: r_model_part.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            r_node.AddDof(VOLUMETRIC_STRAIN, REACTION_STRAIN);
        }

        // Initialize the elements to initialize the constitutive law
        for (auto &r_elem : r_model_part.Elements()) {
            r_elem.Initialize(r_process_info);
        }

        // Types definition
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef MixedGenericCriteria<SparseSpaceType, LocalSpaceType> MixedGenericCriteriaType;

        // Construct the linear solver pointer
        Parameters linear_solver_settings(R"({"solver_type": "skyline_lu_factorization"})");
        LinearSolverFactory<SparseSpaceType, LocalSpaceType> linear_solver_factory;
        auto p_linear_solver = linear_solver_factory.Create(linear_solver_settings);

        // Create the Newton-Raphson strategy
        int max_iterations = 30;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        bool move_mesh_flag = true;
        auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>>();
        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(p_linear_solver);
        Parameters convergence_settings = Parameters(R"({
            "convergence_variables_list" : {
                "pressure" : {
                    "variable"           : "DISPLACEMENT",
                    "relative_tolerance" : 1.0e-3,
                    "absolute_tolerance" : 1.0e-5
                },
                "velocity" : {
                    "variable"           : "VOLUMETRIC_STRAIN",
                    "relative_tolerance" : 1.0e-3,
                    "absolute_tolerance" : 1.0e-5
                }
            }
        })");
        auto p_convergence_criteria = Kratos::make_shared<MixedGenericCriteriaType>(convergence_settings);
        auto p_solving_strategy = Kratos::make_unique<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            r_model_part,
            p_scheme,
            p_linear_solver,
            p_convergence_criteria,
            max_iterations,
            calculate_reactions,
            reform_dof_at_each_iteration,
            move_mesh_flag);
        p_solving_strategy->Check();

        // Set the problem BCs function
        const double tol = 1.0e-8;
        const double disp_x = 1.0;
        const double disp_y = 1.0/4.0;
        array_1d<double,3> aux_disp;
        auto node_bc_func = [&](Node<3>& rNode, double DisplacementFactor){
            // Reset auxiliary displacement vector
            aux_disp = ZeroVector(3);

            // Top boundary
            if (std::abs(rNode.Y0() - 1.0) < tol) {
                aux_disp[1] = -DisplacementFactor * disp_y;
                rNode.Fix(DISPLACEMENT_Y);
                rNode.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            }

            // Bottom boundary
            if (std::abs(rNode.Y0()) < tol) {
                aux_disp[1] = 0.0;
                rNode.Fix(DISPLACEMENT_Y);
                rNode.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            }

            // Left boundary
            if (std::abs(rNode.X0()) < tol) {
                aux_disp[0] = 0.0;
                rNode.Fix(DISPLACEMENT_X);
                rNode.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            }

            // Right boundary
            if (std::abs(rNode.X0() - 1.0) < tol) {
                aux_disp[0] = DisplacementFactor * disp_x;
                rNode.Fix(DISPLACEMENT_X);
                rNode.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            }
        };

        // Set GiD output
        GidIO<> gid_io_str("/media/alm/Data/Works/papers/mixed_disp_J2/TestTotalLagrangianMixedDetJElementBonetPatch", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		gid_io_str.InitializeMesh(0);
		gid_io_str.WriteMesh(r_model_part.GetMesh());
		gid_io_str.FinalizeMesh();
        gid_io_str.InitializeResults(0, r_model_part.GetMesh());

        // Solve the problem incrementally
        const std::size_t n_load_steps = 10000.0;
        for (std::size_t i_step = 1; i_step < n_load_steps + 1; ++i_step) {
            r_model_part.CloneTimeStep(i_step);
            r_model_part.GetProcessInfo().GetValue(STEP) = i_step;
            const double load_factor = double(double(i_step) / double(n_load_steps));
            for (auto& r_node : r_model_part.Nodes()) {
                node_bc_func(r_node, load_factor);
            }
            p_solving_strategy->Solve();
            gid_io_str.WriteNodalResults(DISPLACEMENT, r_model_part.Nodes(), i_step, 0);
            gid_io_str.WriteNodalResults(VOLUMETRIC_STRAIN, r_model_part.Nodes(), i_step, 0);
            gid_io_str.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_VECTOR, r_model_part, i_step);
        }
        gid_io_str.FinalizeResults();

        // // Check results
        // const double tolerance = 1.0e-6;
        // const double expected_vol_strain = 1.4965e-05;
        // const std::vector<double> expected_disp = {-0.000318007, 0.000347937, 0};
        // KRATOS_CHECK_VECTOR_NEAR(p_node_3->FastGetSolutionStepValue(DISPLACEMENT), expected_disp, tolerance);
        // KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN), expected_vol_strain, tolerance);

    }

} // namespace Testing
} // namespace Kratos.
