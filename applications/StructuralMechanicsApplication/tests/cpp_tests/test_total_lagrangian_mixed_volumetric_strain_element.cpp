// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "factories/linear_solver_factory.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
// #include "includes/gid_io.h"
#include "processes/structured_mesh_generator_process.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/total_lagrangian_mixed_volumetric_strain_element.h"

namespace Kratos::Testing
{

    /**
    * Checks the Total Lagrangian mixed Jacobian determinant 2D3N element RHS and LHS
    */
    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedVolumetricStrainElement2D3N, KratosStructuralMechanicsFastSuite)
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
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("TotalLagrangianMixedVolumetricStrainElement2D3N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        // Check RHS and LHS results
        const double tolerance = 1.0e-7;
        const std::vector<double> expected_RHS({58231.8681319,60197.4358974,-0.00451994047619,-19272.1611722,-40923.8095238,0.00872232142857,-38959.7069597,-19273.6263736,0.00227380952381});
        const std::vector<double> expected_LHS_row_0({974614.212454,158304.761905,-325989.010989,-431198.534799,-525956.043956,-325989.010989,-543415.677656,367651.282051,-325989.010989});
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    /**
    * Checks the Total Lagrangian mixed Jacobian determinant 3D4N element RHS and LHS
    */
    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedVolumetricStrainElement3D4N, KratosStructuralMechanicsFastSuite)
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
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 0.0 , 1.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("TotalLagrangianMixedVolumetricStrainElement3D4N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[2] = 0.2;
        p_node_4->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;
        p_node_4->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.05;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(16);
        Matrix LHS = ZeroMatrix(16,16);

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        // Check RHS and LHS results
        const double tolerance = 1.0e-7;
        const std::vector<double> expected_RHS({-9223.89014388,32765.8558739,61472.2259113,-0.0199042607315,9223.89014388,0,0,-0.0121052631599,0,-15913.1975817,-14346.5609408,-0.0104141975643,0,-16852.6582922,-47125.6649705,-0.00559382240399});
        const std::vector<double> expected_LHS_row_0({397802.448164,89460.058478,128948.283247,-63772.4886945,-108460.750823,-130555.762838,-130555.762838,-63772.4886945,-131694.761321,57415.1747145,-19744.1123847,-63732.5294656,-157646.936019,-16319.4703547,21351.591975,-63613.9773907});
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    /**
     * This test implements the simple patch test example in Bonet and Wood book
     */
    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangianMixedVolumetricStrainElementBonetPatch, KratosStructuralMechanicsFastSuite)
    {
        // Skip the test if the constitutive law is not available (i.e. the ConstitutiveLawsApplication is not compiled)
        const std::string claw_name = "HyperElasticPlaneStrain2DLaw";
        if (!KratosComponents<ConstitutiveLaw>::Has(claw_name)) {
            return;
        }

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
		auto p_point_1 = Kratos::make_intrusive<Node>(1, min_x, min_y, 0.0);
		auto p_point_2 = Kratos::make_intrusive<Node>(2, min_x, max_y, 0.0);
		auto p_point_3 = Kratos::make_intrusive<Node>(3, max_x, max_y, 0.0);
		auto p_point_4 = Kratos::make_intrusive<Node>(4, max_x, min_y, 0.0);
		Quadrilateral2D4<Node> square_geometry(p_point_1, p_point_2, p_point_3, p_point_4);
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions" : 2,
            "create_skin_sub_model_part" : false,
			"element_name" : "TotalLagrangianMixedVolumetricStrainElement2D3N"
		})");
		StructuredMeshGeneratorProcess(square_geometry, r_model_part, mesher_parameters).Execute();

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(1);
        p_elem_prop->SetValue(YOUNG_MODULUS, 250.0);
        p_elem_prop->SetValue(POISSON_RATIO, 0.25);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get(claw_name);
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
        bool move_mesh_flag = true;
        bool calculate_reactions = false;
        bool reform_dofs_at_each_iteration = false;
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
            reform_dofs_at_each_iteration,
            move_mesh_flag);
        p_solving_strategy->Check();
        p_solving_strategy->SetEchoLevel(0);
        p_convergence_criteria->SetEchoLevel(0);

        // Set the problem BCs function
        const double tol = 1.0e-8;
        const double disp_x = 1.0;
        const double disp_y = 1.0/4.0;
        array_1d<double,3> aux_disp;
        auto node_bc_func = [&](Node& rNode, double DisplacementFactor){
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

        // // Set GiD output
        // GidIO<> gid_io_str("/media/alm/Data/Works/papers/mixed_disp_J2/TestTotalLagrangianMixedVolumetricStrainElementBonetPatch", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		// gid_io_str.InitializeMesh(0);
		// gid_io_str.WriteMesh(r_model_part.GetMesh());
		// gid_io_str.FinalizeMesh();
        // gid_io_str.InitializeResults(0, r_model_part.GetMesh());

        // Solve the problem incrementally
        const std::size_t n_load_steps = 1.0;
        for (std::size_t i_step = 1; i_step < n_load_steps + 1; ++i_step) {
            r_model_part.CloneTimeStep(i_step);
            r_model_part.GetProcessInfo().GetValue(STEP) = i_step;
            const double load_factor = double(double(i_step) / double(n_load_steps));
            for (auto& r_node : r_model_part.Nodes()) {
                node_bc_func(r_node, load_factor);
            }
            p_solving_strategy->Solve();
            // gid_io_str.WriteNodalResults(DISPLACEMENT, r_model_part.Nodes(), i_step, 0);
            // gid_io_str.WriteNodalResults(VOLUMETRIC_STRAIN, r_model_part.Nodes(), i_step, 0);
            // gid_io_str.PrintOnGaussPoints(PK2_STRESS_VECTOR, r_model_part, i_step);
            // gid_io_str.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_VECTOR, r_model_part, i_step);
        }
        // gid_io_str.FinalizeResults();

        // Check results
        const double tolerance = 1.0e-6;
        const double expected_vol_strain = 0.5;
        const std::vector<double> expected_E = {1.5,-0.218750000047,0.0};
        const std::vector<double> expected_PK2 = {85.1366277052,-5.69509189368,0.0};
        for (const auto& r_node : r_model_part.Nodes()) {
            KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN), expected_vol_strain, tolerance);
        }
        for (auto& r_elem : r_model_part.Elements()) {
            std::vector<Vector> obtained_E;
            std::vector<Vector> obtained_PK2;
            r_elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, obtained_PK2, r_model_part.GetProcessInfo());
            r_elem.CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, obtained_E, r_model_part.GetProcessInfo());
            for (std::size_t i_gauss = 0; i_gauss < obtained_E.size(); ++i_gauss) {
                KRATOS_EXPECT_VECTOR_NEAR(obtained_E[i_gauss], expected_E, tolerance);
                KRATOS_EXPECT_VECTOR_NEAR(obtained_PK2[i_gauss], expected_PK2, tolerance);
            }
        }
    }

} // namespace Kratos::Testing
