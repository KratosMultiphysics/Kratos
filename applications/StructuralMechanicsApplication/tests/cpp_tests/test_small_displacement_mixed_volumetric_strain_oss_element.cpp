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
// #include "includes/gid_io.h"
#include "includes/variables.h"
#include "factories/linear_solver_factory.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "structural_mechanics_fast_suite.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_element.h"

namespace Kratos::Testing
{

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainOssElement2D3N, KratosStructuralMechanicsFastSuite)
{

    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

    auto& r_process_info = r_model_part.GetProcessInfo();

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);
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
    auto p_node_3 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
    std::vector<ModelPart::IndexType> element_nodes {1,2,3};
    auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 1, element_nodes, p_elem_prop);

    // Set a fake displacement and volumetric strain field to compute the residual
    array_1d<double, 3> aux_disp = ZeroVector(3);
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    aux_disp[1] = 0.1;
    p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
    p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
    p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;

    // Set fake projection values to calculate the residual
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_3->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 1e-3;
    p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 2e-3;
    p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 3e-3;

    // Compute RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    // Check RHS and LHS results
    const double tolerance = 1.0e-5;
    const std::vector<double> expected_RHS({51173.0769231, 51173.0769231, -1816.54062688, -12711.5384615, -38461.5384615, 10555.3527785, -38461.5384615, -12711.5384615, 3972.72630992});
    const std::vector<double> expected_LHS_row_0({778846.153846, 9615.38461538, -317307.692308, -394230.769231, -384615.384615, -317307.692308, -384615.384615, 375000, -317307.692308});
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainOssElement2D3NDynamic, KratosStructuralMechanicsFastSuite)
{

    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DELTA_TIME] = 0.1;

    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_PROJECTION);
    r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN_PROJECTION);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e06);
    p_elem_prop->SetValue(POISSON_RATIO, 0.3);
    p_elem_prop->SetValue(DENSITY, 1.0e3);
    const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
    std::vector<ModelPart::IndexType> element_nodes {1,2,3};
    auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 1, element_nodes, p_elem_prop);

    // Set a fake displacement, volumetric strain and acceleration fields to compute the residual
    array_1d<double, 3> aux_disp = ZeroVector(3);
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
    aux_disp[1] = 0.1;
    p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;

    array_1d<double,3> aux_acc = ZeroVector(3);
    p_node_1->FastGetSolutionStepValue(ACCELERATION) = aux_acc;
    aux_acc[0] = 100.0;
    p_node_2->FastGetSolutionStepValue(ACCELERATION) = aux_acc;
    aux_acc[1] = 200.0;
    p_node_3->FastGetSolutionStepValue(ACCELERATION) = aux_acc;

    p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
    p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
    p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;

    // Set fake projection values to calculate the residual
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_3->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION) = aux_disp;
    p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 1e-3;
    p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 2e-3;
    p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION) = 3e-3;

    // Compute RHS, LHS and mass matrix
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);
    Matrix MassMatrix = ZeroMatrix(9,9);

    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    p_element->FinalizeSolutionStep(r_process_info); // Fake call to the FinalizeSolutionStep in order to upgrade the subscales dynamic component
    p_element->FinalizeSolutionStep(r_process_info); // Fake call to the FinalizeSolutionStep in order to upgrade the subscales dynamic component
    p_element->CalculateLocalSystem(LHS, RHS, r_process_info);
    p_element->CalculateMassMatrix(MassMatrix, r_process_info);

    // Check RHS and LHS results
    const double tolerance = 1.0e-5;
    const std::vector<double> expected_RHS({51181.5777242, 51173.076373, -1831.8931645, -12703.0376604, -38461.5390116, 10570.705396, -38453.0376604, -12711.5390116, 3972.72623009});
    const std::vector<double> expected_LHS_row_0({778846.153846, 9615.38461538, -327036.58062, -394230.769231, -384615.384615, -307578.803995, -384615.384615, 375000, -317307.692308});
    const std::vector<double> expected_mass_row_0({83.3333333333, 0, 0, 41.6666666667, 0, 0, 41.6666666667, 0, 0});
    const std::vector<double> expected_mass_row_2({0, 0, 0, 0, 0, 0, 0, 0, 0});
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(row(MassMatrix,0), expected_mass_row_0, tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(row(MassMatrix,2), expected_mass_row_2, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedVolumetricStrainOssElementZienkiewiczPatch, KratosStructuralMechanicsFastSuite)
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

    // Activate the OSS projections
    r_process_info.GetValue(OSS_SWITCH) = true;

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
    auto p_element_1 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 1, element_nodes_1, p_elem_prop);
    auto p_element_2 = r_model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 2, element_nodes_2, p_elem_prop);

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
    p_solving_strategy->SetEchoLevel(0);
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
    const double expected_vol_strain = 1.49650698603e-05;
    const double expected_vol_strain_proj = -0.000531463341208;
    const std::vector<double> expected_disp = {-0.000533062, 0.000562992, 0.0};
    const std::vector<double> expected_disp_proj = {-2.20921422407, 2.20921422407, 0.0};
    KRATOS_CHECK_VECTOR_NEAR(p_node_3->FastGetSolutionStepValue(DISPLACEMENT), expected_disp, tolerance);
    KRATOS_CHECK_VECTOR_NEAR(p_node_3->FastGetSolutionStepValue(DISPLACEMENT_PROJECTION), expected_disp_proj, tolerance);
    KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN), expected_vol_strain, tolerance);
    KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION), 0.0, tolerance);
    KRATOS_CHECK_NEAR(p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION), expected_vol_strain_proj, tolerance);

    // // GiD output
    // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/test_small_displacement_mixed_volumetric_strain_oss_element_zienkiewicz_patch", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_fluid.InitializeMesh(0.0);
    // gid_io_fluid.WriteMesh(r_model_part.GetMesh());
    // gid_io_fluid.FinalizeMesh();
    // gid_io_fluid.InitializeResults(0, r_model_part.GetMesh());
    // gid_io_fluid.WriteNodalResults(DISPLACEMENT, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.WriteNodalResults(VOLUMETRIC_STRAIN, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.WriteNodalResults(DISPLACEMENT_PROJECTION, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.WriteNodalResults(VOLUMETRIC_STRAIN_PROJECTION, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.FinalizeResults();

}

} // namespace Kratos::Testing
