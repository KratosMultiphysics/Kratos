//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "includes/gid_io.h"
#include "processes/structured_mesh_generator_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/two_fluid_history_projection_utility.h"

namespace Kratos {

namespace Testing {

void SetModelPart(ModelPart& rNewModelPart)
{

    rNewModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rNewModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rNewModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Generate nodes and elements
    Properties::Pointer p_properties(new Properties(0));
    auto p_point_1 = rNewModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_point_2 = rNewModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_point_3 = rNewModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
    auto p_point_4 = rNewModelPart.CreateNewNode(4, 3.0, 0.0, 0.0);
    auto p_point_5 = rNewModelPart.CreateNewNode(5, 4.0, 0.0, 0.0);
    auto p_point_6 = rNewModelPart.CreateNewNode(6, 0.0, 1.0, 0.0);
    auto p_point_7 = rNewModelPart.CreateNewNode(7, 1.0, 1.0, 0.0);
    auto p_point_8 = rNewModelPart.CreateNewNode(8, 2.0, 1.0, 0.0);
    auto p_point_9 = rNewModelPart.CreateNewNode(9, 3.0, 1.0, 0.0);
    auto p_point_10 = rNewModelPart.CreateNewNode(10, 4.0, 1.0, 0.0);
    rNewModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 6}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 2, {2, 7, 6}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 3, {2, 3, 7}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 4, {3, 8, 7}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 5, {3, 4, 8}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 6, {4, 9, 8}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 7, {4, 5, 9}, p_properties);
    rNewModelPart.CreateNewElement("Element2D3N", 8, {5, 10, 9}, p_properties);
}

/**
 * Checks the embedded skin visualization process for a unique triangle with the standard shape functions
 */
KRATOS_TEST_CASE_IN_SUITE(TwoFluidHistoryProjectionConstantVelocity, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model model;
    ModelPart& main_model_part = model.CreateModelPart("MainModelPart");
    SetModelPart(main_model_part);
    main_model_part.SetBufferSize(2);
    main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    main_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.1);

    // Set distance and velocity fields
    array_1d<double,3> aux_v({1.0,0.0,0.0});
    for (auto& node : main_model_part.Nodes()) {
        node.FastGetSolutionStepValue(DISTANCE, 0) = node.X() - 2.5;
        node.FastGetSolutionStepValue(DISTANCE, 1) = node.X() - 1.5;
        node.FastGetSolutionStepValue(VELOCITY, 0) = aux_v;
        node.FastGetSolutionStepValue(VELOCITY, 1) = aux_v;
    }

    // Call the history projection utility
    const double search_factor = 4.0;
    const bool compute_nodal_h = true;
    const double particle_layer_thickness = 1.5;
    for (auto &node : main_model_part.Nodes())
    {
        KRATOS_WATCH(node.FastGetSolutionStepValue(DISTANCE, 0));
        KRATOS_WATCH(node.FastGetSolutionStepValue(DISTANCE, 1));
    }
    TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        main_model_part,
        compute_nodal_h,
        particle_layer_thickness,
        search_factor);

    GidIO<> gid_io_fluid("/home/uchasco/Desktop/test_two_fluid_history_projection_utility", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_fluid.InitializeMesh(0.0);
    gid_io_fluid.WriteMesh(main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(0, main_model_part.GetMesh());
    gid_io_fluid.WriteNodalResults(DISTANCE, main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, main_model_part.Nodes(), 0);

    gid_io_fluid.FinalizeResults();

    // Check values
    const double tolerance = 1.0e-8;
    for (auto& r_node : main_model_part.Nodes()) {
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(VELOCITY, 0))
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(VELOCITY, 1))
        auto& r_v_mesh = r_node.FastGetSolutionStepValue(MESH_VELOCITY);
        std::cout << "Node: " << r_node.Id() << " v_mesh: [" << r_v_mesh[0] << "," <<  r_v_mesh[1] << "," << r_v_mesh[2] << "]" << std::endl;
        KRATOS_CHECK_VECTOR_NEAR(r_node.FastGetSolutionStepValue(VELOCITY, 0), aux_v, tolerance);
        KRATOS_CHECK_VECTOR_NEAR(r_node.FastGetSolutionStepValue(VELOCITY, 1), aux_v, tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TwoFluidHistoryProjectionSquareConstantVelocityLinear, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model model;
    ModelPart &r_main_model_part = model.CreateModelPart("MainModelPart");
    r_main_model_part.SetBufferSize(2);
    r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    r_main_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.1);

    // Add required variables
    r_main_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_main_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_main_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Create the test geometry
    auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions": 14,
        "element_name": "Element2D3N"
    })");
    StructuredMeshGeneratorProcess(geometry, r_main_model_part, mesher_parameters).Execute();

    // Set distance and velocity fields
    for (auto& r_node : r_main_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISTANCE) = r_node.X() - 0.5;
        r_node.FastGetSolutionStepValue(DISTANCE, 1) = r_node.X() - 0.5;
        r_node.FastGetSolutionStepValue(VELOCITY) = array_1d<double,3>({1.0,0.0,0.0});
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = array_1d<double,3>({1.0,0.0,0.0});
    }

    GidIO<> gid_io_fluid("/home/uchasco/Desktop/ProjectionUtilityCppTests/test_two_fluid_history_projection_utility_square_constant_velocity_linear", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_fluid.InitializeMesh(0.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(0, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 0);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 0);
    gid_io_fluid.FinalizeResults();

    // Do the particle projection
    const double search_factor = 10;
    const bool compute_nodal_h = true;
    const double particle_layer_thickness = 0.25;
    TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        r_main_model_part,
        compute_nodal_h,
        particle_layer_thickness,
        search_factor);

    gid_io_fluid.InitializeMesh(1.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(1, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 1);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 1);
    gid_io_fluid.FinalizeResults();
}

KRATOS_TEST_CASE_IN_SUITE(TwoFluidHistoryProjectionSquareConstantVelocityLinearAngle, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model model;
    ModelPart &r_main_model_part = model.CreateModelPart("MainModelPart");
    r_main_model_part.SetBufferSize(2);
    r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    r_main_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.1);

    // Add required variables
    r_main_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_main_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_main_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Create the test geometry
    auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions": 28,
        "element_name": "Element2D3N"
    })");
    StructuredMeshGeneratorProcess(geometry, r_main_model_part, mesher_parameters).Execute();

    // Set distance and velocity fields
    for (auto &r_node : r_main_model_part.Nodes())
    {
        r_node.FastGetSolutionStepValue(DISTANCE) = -r_node.Y()-r_node.X() +1.0;
        r_node.FastGetSolutionStepValue(DISTANCE, 1) = -r_node.Y() - r_node.X() + 1.0;
        r_node.FastGetSolutionStepValue(VELOCITY) = array_1d<double, 3>({1.0, 1.0, 0.0});
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = array_1d<double, 3>({1.0, 1.0, 0.0});
    }

    GidIO<> gid_io_fluid("/home/uchasco/Desktop/ProjectionUtilityCppTests/test_two_fluid_history_projection_utility_square_constant_velocity_linear_angle", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_fluid.InitializeMesh(0.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(0, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 0);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 0);
    gid_io_fluid.FinalizeResults();

    // Do the particle projection
    const double search_factor = 4.0;
    const bool compute_nodal_h = true;
    const double particle_layer_thickness = 20;
    TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        r_main_model_part,
        compute_nodal_h,
        particle_layer_thickness,
        search_factor);

    gid_io_fluid.InitializeMesh(1.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(1, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 1);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 1);
    gid_io_fluid.FinalizeResults();
}

KRATOS_TEST_CASE_IN_SUITE(TwoFluidHistoryProjectionSquareConstantVelocitySinusoidal, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model model;
    ModelPart &r_main_model_part = model.CreateModelPart("MainModelPart");
    r_main_model_part.SetBufferSize(2);
    r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    r_main_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.1);

    // Add required variables
    r_main_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_main_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_main_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Create the test geometry
    auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions": 7,
        "element_name": "Element2D3N"
    })");
    StructuredMeshGeneratorProcess(geometry, r_main_model_part, mesher_parameters).Execute();

    // Set distance and velocity fields
    for (auto& r_node : r_main_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISTANCE) = r_node.X() - (0.5 + 0.1 * std::cos(Globals::Pi * (1.0 - r_node.Y())));
        r_node.FastGetSolutionStepValue(DISTANCE, 1) = r_node.X() - (0.5 + 0.1 * std::cos(Globals::Pi * (1.0 - r_node.Y())));
        r_node.FastGetSolutionStepValue(VELOCITY) = array_1d<double,3>({1.0,0.0,0.0});
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = array_1d<double,3>({1.0,0.0,0.0});
    }

    GidIO<> gid_io_fluid("/home/uchasco/Desktop/ProjectionUtilityCppTests/test_two_fluid_history_projection_utility_square_constant_velocity_sinusoidal", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_fluid.InitializeMesh(0.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(0, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 0);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 0);
    gid_io_fluid.FinalizeResults();

    // Do the particle projection
    const double search_factor = 4.0;
    const bool compute_nodal_h = true;
    const double particle_layer_thickness = 0.25;
    TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        r_main_model_part,
        compute_nodal_h,
        particle_layer_thickness,
        search_factor);

    gid_io_fluid.InitializeMesh(1.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(1, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 1);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 1);
    gid_io_fluid.FinalizeResults();
}

KRATOS_TEST_CASE_IN_SUITE(TwoFluidHistoryProjectionSquareConstantVelocityCircle, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model model;
    ModelPart &r_main_model_part = model.CreateModelPart("MainModelPart");
    r_main_model_part.SetBufferSize(2);
    r_main_model_part.GetProcessInfo().SetValue(TIME, 0.0);
    r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    r_main_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.3333333333333333);

    // Add required variables
    r_main_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_main_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_main_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Create the test geometry
    auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions": 15,
        "element_name": "Element2D3N"
    })");
    StructuredMeshGeneratorProcess(geometry, r_main_model_part, mesher_parameters).Execute();

    // Set distance and velocity fields
    for (auto& r_node : r_main_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISTANCE) = std::sqrt(std::pow(r_node.X() - 0.25, 2) + std::pow(r_node.Y() - 0.25, 2)) - 0.125;
        r_node.FastGetSolutionStepValue(DISTANCE,1) = std::sqrt(std::pow(r_node.X() - 0.25, 2) + std::pow(r_node.Y() - 0.25, 2)) - 0.125;
        r_node.FastGetSolutionStepValue(VELOCITY) = array_1d<double,3>({1.0,1.0,0.0});
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = array_1d<double,3>({1.0,1.0,0.0});
    }

    GidIO<> gid_io_fluid("/home/uchasco/Desktop/ProjectionUtilityCppTests/test_two_fluid_history_projection_utility_square_constant_velocity_circle", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_fluid.InitializeMesh(0.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(0, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 0);
    gid_io_fluid.WriteNodalFlags(FREE_SURFACE, "FREE_SURFACE", r_main_model_part.Nodes(), 0);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 0, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 0);
    gid_io_fluid.FinalizeResults();

    // Do the particle projection
    const std::size_t steps =1 ;
    const double dt = r_main_model_part.GetProcessInfo().GetValue(DELTA_TIME);
    r_main_model_part.GetProcessInfo().SetValue(DELTA_TIME, dt / steps);

    for (std::size_t i = 0; i < steps; ++i) {
        KRATOS_WATCH(i)
        const double time = r_main_model_part.GetProcessInfo().GetValue(TIME);
        const double dt = r_main_model_part.GetProcessInfo().GetValue(DELTA_TIME);
        r_main_model_part.CloneTimeStep(time + dt);
        const double search_factor = 4;
        const bool compute_nodal_h = true;
        const double particle_layer_thickness = 0.2;
        TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
            r_main_model_part,
            compute_nodal_h,
            particle_layer_thickness,
            search_factor);
    }

    gid_io_fluid.InitializeMesh(1.0);
    gid_io_fluid.WriteMesh(r_main_model_part.GetMesh());
    gid_io_fluid.FinalizeMesh();
    gid_io_fluid.InitializeResults(1, r_main_model_part.GetMesh());
    gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", r_main_model_part.Nodes(), 1);
    gid_io_fluid.WriteNodalFlags(FREE_SURFACE, "FREE_SURFACE", r_main_model_part.Nodes(), 1);
    gid_io_fluid.WriteNodalResults(DISTANCE, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResults(MESH_VELOCITY, r_main_model_part.Nodes(), 1, 0);
    gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_main_model_part.Nodes(), 1);
    gid_io_fluid.FinalizeResults();
}

} // namespace Testing

}  // namespace Kratos.
