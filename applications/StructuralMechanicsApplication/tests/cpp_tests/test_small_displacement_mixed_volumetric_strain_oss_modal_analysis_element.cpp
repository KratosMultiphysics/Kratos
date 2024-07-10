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
#include "structural_mechanics_fast_suite.h"
#include "includes/gid_io.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_element.h"
#include "factories/linear_solver_factory.h"
#include "includes/cfd_variables.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos::Testing
{

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainOssModalAnalysisElement2D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_PROJECTION);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN_PROJECTION);

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
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainModalAnalysisElement2D4N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp = ZeroVector(3);
        noalias(p_node_1->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        noalias(p_node_2->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        aux_disp[1] = 0.1;
        noalias(p_node_3->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        noalias(p_node_4->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;
        p_node_4->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.03;

        // Set the body force
        array_1d<double,3> body_force = ZeroVector(3);
        noalias(p_node_1->FastGetSolutionStepValue(VOLUME_ACCELERATION)) = body_force;
        noalias(p_node_2->FastGetSolutionStepValue(VOLUME_ACCELERATION)) = body_force;
        noalias(p_node_3->FastGetSolutionStepValue(VOLUME_ACCELERATION)) = body_force;
        noalias(p_node_4->FastGetSolutionStepValue(VOLUME_ACCELERATION)) = body_force;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(24);
        Matrix LHS = ZeroMatrix(24,24);

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({69368.1, 145833, -26251.5, -69368.1, 146062, -34188, -68681.3, -146062, 21596.5, 68681.3, -145833, 61507.9, -0.00208333, 0.0104167, 46398, -0.00208333, 0.00833333, 47008.5, -0.00416667, 0.00833333, 44566.5, -0.00416667, 0.0104167, 43345.5});
        const std::vector<double> expected_LHS_row_0({1.24542e+06, 549451, 45787.5, -860806, 164835, 45787.5, -622711, -549451, 22893.8, 238095, -164835, 22893.8, 0, 0, 366300, 0, 0, 366300, 0, 0, 183150, 0, 0, 183150});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * Zienkiewicz patch test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainOssModalAnalysisElementZienkiewiczPatch, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart",1);
        auto& r_process_info = r_model_part.GetProcessInfo();
        r_process_info.GetValue(DOMAIN_SIZE) = 2;

        r_model_part.AddNodalSolutionStepVariable(REACTION);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_PROJECTION);
        r_model_part.AddNodalSolutionStepVariable(REACTION_STRAIN);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN_PROJECTION);

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
        auto p_element_1 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainModalAnalysisElement2D3N", 1, element_nodes_1, p_elem_prop);
        auto p_element_2 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainModalAnalysisElement2D3N", 2, element_nodes_2, p_elem_prop);

        // Create the load condition
        auto p_cond_prop = r_model_part.CreateNewProperties(0);
        std::vector<ModelPart::IndexType> condition_nodes_1 {3};
        auto p_condition_1 = r_model_part.CreateNewCondition("PointLoadCondition2D1N", 1, condition_nodes_1, p_cond_prop);
        array_1d<double,3> point_load = ZeroVector(3);
        point_load[1] = 1.0;
        p_condition_1->SetValue(POINT_LOAD, point_load);

        // Add DOFs
        for (auto &r_node: r_model_part.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
            r_node.AddDof(VOLUMETRIC_STRAIN);
            r_node.AddDof(DISPLACEMENT_PROJECTION_X);
            r_node.AddDof(DISPLACEMENT_PROJECTION_Y);
            r_node.AddDof(DISPLACEMENT_PROJECTION_Z);
            r_node.AddDof(VOLUMETRIC_STRAIN_PROJECTION);
        }

        // Initialize the elements to initialize the constitutive law
        for (auto &r_elem : r_model_part.Elements()) {
            r_elem.Initialize(r_process_info);
        }

        // Construct the linear solver pointer
        Parameters linear_solver_settings(R"({"solver_type": "skyline_lu_factorization"})");
        LinearSolverFactory<SparseSpaceType, LocalSpaceType> linear_solver_factory;
        auto p_linear_solver = linear_solver_factory.Create(linear_solver_settings);

        // Create the non-linear strategy
        Scheme<SparseSpaceType,LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>>();

        const std::size_t max_it = 20;
        const bool move_mesh_flag = true;
        const bool calculate_reactions = false;
        const bool reform_dof_at_each_iteration = false;
        auto p_conv_criteria = Kratos::make_shared<DisplacementCriteria<SparseSpaceType, LocalSpaceType>>(1e-10, 1e-9);
        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(p_linear_solver);
        auto p_solving_strategy = Kratos::make_unique<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            r_model_part,
            p_scheme,
            p_conv_criteria,
            p_builder_and_solver,
            max_it,
            calculate_reactions,
            reform_dof_at_each_iteration,
            move_mesh_flag);
        // p_solving_strategy->SetEchoLevel(0);
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

        // // Check results
        // const double tolerance = 1.0e-6;
        // const double expected_vol_strain = 1.4965e-05;
        // const std::vector<double> expected_disp = {-0.000318007, 0.000347937, 0};
        // KRATOS_CHECK_VECTOR_NEAR(p_node_3->FastGetSolutionStepValue(DISPLACEMENT), expected_disp, tolerance);
        // KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN), expected_vol_strain, tolerance);

        // GiD output
        GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/test_small_displacement_mixed_volumetric_strain_oss_modal_analysis_element_zienkiewicz_patch", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		gid_io_fluid.InitializeMesh(0.0);
		gid_io_fluid.WriteMesh(r_model_part.GetMesh());
		gid_io_fluid.FinalizeMesh();
		gid_io_fluid.InitializeResults(0, r_model_part.GetMesh());
		gid_io_fluid.WriteNodalResults(DISPLACEMENT, r_model_part.Nodes(), 0, 0);
		gid_io_fluid.WriteNodalResults(VOLUMETRIC_STRAIN, r_model_part.Nodes(), 0, 0);
		gid_io_fluid.WriteNodalResults(DISPLACEMENT_PROJECTION, r_model_part.Nodes(), 0, 0);
		gid_io_fluid.WriteNodalResults(VOLUMETRIC_STRAIN_PROJECTION, r_model_part.Nodes(), 0, 0);
		gid_io_fluid.FinalizeResults();

    }

} // namespace Kratos::Testing
