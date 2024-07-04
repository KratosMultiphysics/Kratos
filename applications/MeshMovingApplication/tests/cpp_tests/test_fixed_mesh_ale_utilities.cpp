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

// Project includes
#include "containers/model.h"
#include "includes/expect.h"
// #include "includes/gid_io.h"
#include "includes/mesh_moving_variables.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/structured_mesh_generator_process.h"

// Application includes
#include "custom_utilities/fixed_mesh_ale_utilities.h"
#include "tests/cpp_tests/mesh_moving_fast_suite.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(FixedMeshALEUtilities2D, MeshMovingApplicationFastSuite)
{
    Model current_model;

    // Generate the origin model part (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node > geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 7,
        "element_name": "Element2D3N"
    })");

    ModelPart& origin_model_part = current_model.CreateModelPart("OriginModelPart");
    origin_model_part.SetBufferSize(3);
    origin_model_part.AddNodalSolutionStepVariable(DISTANCE);
    origin_model_part.AddNodalSolutionStepVariable(VELOCITY);
    origin_model_part.AddNodalSolutionStepVariable(PRESSURE);
    origin_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    origin_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    origin_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
    StructuredMeshGeneratorProcess(geometry, origin_model_part, mesher_parameters).Execute();

    // Fix the boundary mesh displacement
    const double max_x = 1.0;
    const double min_x = 0.0;
    const double max_y = 1.0;
    const double min_y = 0.0;
    const double coord_tol = 1e-3;
    for (auto it_node : origin_model_part.NodesArray()){
        if ((std::abs(it_node->X() - min_x) < coord_tol) || (std::abs(it_node->X() - max_x) < coord_tol)) {
            it_node->Fix(MESH_DISPLACEMENT_X);
        }
        if ((std::abs(it_node->Y() - min_y) < coord_tol) || (std::abs(it_node->Y() - max_y) < coord_tol)) {
            it_node->Fix(MESH_DISPLACEMENT_Y);
        }
    }

    // Create a fake time loop to fill the buffer
    const double n_steps = 3;
    const double delta_time = 0.1;
    for (unsigned int i_step = 0; i_step < n_steps; ++i_step){
        origin_model_part.CloneTimeStep(i_step * delta_time);
        double p_val = i_step * delta_time;
        array_1d<double,3> v_val = ZeroVector(3);
        for (auto it_node : origin_model_part.NodesArray()){
            it_node->GetSolutionStepValue(PRESSURE) = p_val;
            v_val(0) = i_step * it_node->X();
            it_node->GetSolutionStepValue(VELOCITY) = v_val;
            it_node->GetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        }
    }

    // Set the virtual model part
    ModelPart& virtual_model_part = current_model.CreateModelPart("VirtualModelPart");

    // Set the structure model part
    ModelPart& str_model_part =current_model.CreateModelPart("StructureModelPart");
    str_model_part.SetBufferSize(3);
    str_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);
    str_model_part.CreateNewNode(1000, 0.6, 0.2, 0.0);
    str_model_part.CreateNewNode(2000, 0.6, 0.4, 0.0);
    str_model_part.CreateNewNode(3000, 0.8, 0.4, 0.0);
    str_model_part.CreateNewNode(4000, 0.8, 0.2, 0.0);
    str_model_part.CreateNewElement("Element2D2N", 1, {{1000,2000}}, p_prop);
    str_model_part.CreateNewElement("Element2D2N", 2, {{2000,3000}}, p_prop);
    str_model_part.CreateNewElement("Element2D2N", 3, {{3000,4000}}, p_prop);
    str_model_part.CreateNewElement("Element2D2N", 4, {{4000,1000}}, p_prop);

    // Set the structure mesh movement
    array_1d<double,3> str_mov = ZeroVector(3);
    for (unsigned int i_buffer = 0; i_buffer < 3; ++i_buffer) {
        for (auto it_str_node : str_model_part.NodesArray()){
            str_mov(0) = -(it_str_node->X()) * 0.05 * (3 - i_buffer);
            str_mov(1) = (it_str_node->Y()) * 0.1 * (3 - i_buffer);
            it_str_node->FastGetSolutionStepValue(DISPLACEMENT, i_buffer) = str_mov;
        }
    }

    // Compute distance
    CalculateDistanceToSkinProcess<2>(origin_model_part, str_model_part).Execute();

    // Set the FM-ALE utility
    Parameters fm_ale_settings(R"({
        "virtual_model_part_name" : "VirtualModelPart",
        "structure_model_part_name" : "StructureModelPart"
    })");
    auto p_mesh_moving = Kratos::make_shared<FixedMeshALEUtilities>(current_model, fm_ale_settings);

    // Fill the virtual model part geometry
    p_mesh_moving->Initialize(origin_model_part);

    // Fill the virtual mesh values
    p_mesh_moving->SetVirtualMeshValuesFromOriginMesh();

    // Execute the FM-ALE operations
    const unsigned int buffer_size = 3;
    p_mesh_moving->ComputeMeshMovement(delta_time);
    p_mesh_moving->ProjectVirtualValues<2>(origin_model_part, buffer_size);
    p_mesh_moving->UndoMeshMovement();

    // GidIO<> gid_io_structure("/home/rzorrilla/Desktop/structure_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_structure.InitializeMesh(0.00);
    // gid_io_structure.WriteMesh(str_model_part.GetMesh());
    // gid_io_structure.FinalizeMesh();
    // gid_io_structure.InitializeResults(0, str_model_part.GetMesh());
    // gid_io_structure.WriteNodalResults(DISPLACEMENT, str_model_part.Nodes(), 0, 0);
    // gid_io_structure.FinalizeResults();

    // GidIO<> gid_io_origin("/home/rzorrilla/Desktop/origin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_origin.InitializeMesh(0.00);
    // gid_io_origin.WriteMesh(origin_model_part.GetMesh());
    // gid_io_origin.FinalizeMesh();
    // gid_io_origin.InitializeResults(0, origin_model_part.GetMesh());
    // gid_io_origin.WriteNodalResults(DISTANCE, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.WriteNodalResults(VELOCITY, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.WriteNodalResults(PRESSURE, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.WriteNodalResults(DISPLACEMENT, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.WriteNodalResults(MESH_VELOCITY, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.WriteNodalResults(MESH_DISPLACEMENT, origin_model_part.Nodes(), 0, 0);
    // gid_io_origin.FinalizeResults();

    // GidIO<> gid_io_virtual("/home/rzorrilla/Desktop/virtual_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_virtual.InitializeMesh(0.00);
    // gid_io_virtual.WriteMesh(virtual_model_part.GetMesh());
    // gid_io_virtual.FinalizeMesh();
    // gid_io_virtual.InitializeResults(0, virtual_model_part.GetMesh());
    // gid_io_virtual.WriteNodalResults(VELOCITY, virtual_model_part.Nodes(), 0, 0);
    // gid_io_virtual.WriteNodalResults(PRESSURE, virtual_model_part.Nodes(), 0, 0);
    // gid_io_virtual.WriteNodalResults(DISPLACEMENT, virtual_model_part.Nodes(), 0, 0);
    // gid_io_virtual.WriteNodalResults(MESH_VELOCITY, virtual_model_part.Nodes(), 0, 0);
    // gid_io_virtual.WriteNodalResults(MESH_DISPLACEMENT, virtual_model_part.Nodes(), 0, 0);
    // gid_io_virtual.FinalizeResults();

    // Check the obtained results
    const double tol = 1e-5;

    // Check that the intersected elements nodes respect the fixity
    const auto node_orig_51_u_mesh = origin_model_part.pGetNode(64)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
    const auto node_virt_51_u_mesh = virtual_model_part.pGetNode(64)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
    for (std::size_t i = 0; i < 3; ++i) {
        KRATOS_EXPECT_NEAR(node_orig_51_u_mesh[i], node_virt_51_u_mesh[i], tol);
    }

    // Check the obtained displacement values in the virtual mesh
    const auto u_mesh_29 = virtual_model_part.pGetNode(29)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
    const auto u_mesh_53 = virtual_model_part.pGetNode(54)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
    const std::array<double,6> expected_values_u_mesh{{-0.0292061,0.0299231,0,-0.0163507,0.0340681,0}};
    const std::array<double,6> obtained_values_u_mesh{{u_mesh_29[0], u_mesh_29[1], u_mesh_29[2], u_mesh_53[0], u_mesh_53[1], u_mesh_53[2]}};
    for (std::size_t i = 0; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(obtained_values_u_mesh[i], expected_values_u_mesh[i], tol);
    }

    // Check the projected mesh velocity in the origin mesh
    const auto v_mesh_29 = origin_model_part.pGetNode(29)->FastGetSolutionStepValue(MESH_VELOCITY);
    const auto v_mesh_53 = origin_model_part.pGetNode(54)->FastGetSolutionStepValue(MESH_VELOCITY);
    const std::array<double,6> expected_values_v_mesh{{-0.30534,0.351295,0,-0.282042,0.311169,0}};
    const std::array<double,6> obtained_values_v_mesh{{v_mesh_29[0], v_mesh_29[1], v_mesh_29[2], v_mesh_53[0], v_mesh_53[1], v_mesh_53[2]}};
    for (std::size_t i = 0; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(obtained_values_v_mesh[i], expected_values_v_mesh[i], tol);
    }

    // Check the projected values in the origin mesh
    const auto v_29 = origin_model_part.pGetNode(29)->FastGetSolutionStepValue(VELOCITY);
    const auto v_53 = origin_model_part.pGetNode(54)->FastGetSolutionStepValue(VELOCITY);
    const auto p_29 = origin_model_part.pGetNode(29)->FastGetSolutionStepValue(PRESSURE);
    const auto p_53 = origin_model_part.pGetNode(54)->FastGetSolutionStepValue(PRESSURE);
    const std::array<double,8> expected_projected_values{{0.857143,0,0, 0.2,1.71429,0,0,0.2}};
    const std::array<double,8> obtained_projected_values{{v_29[0], v_29[1], v_29[2], p_29, v_53[0], v_53[1], v_53[2], p_53}};
    for (std::size_t i = 0; i < 8; ++i) {
        KRATOS_EXPECT_NEAR(obtained_projected_values[i], expected_projected_values[i], tol);
    }

    const auto v_n_29 = origin_model_part.pGetNode(29)->FastGetSolutionStepValue(VELOCITY, 1);
    const auto v_n_53 = origin_model_part.pGetNode(54)->FastGetSolutionStepValue(VELOCITY, 1);
    const auto p_n_29 = origin_model_part.pGetNode(29)->FastGetSolutionStepValue(PRESSURE, 1);
    const auto p_n_53 = origin_model_part.pGetNode(54)->FastGetSolutionStepValue(PRESSURE, 1);
    const std::array<double,8> expected_projected_values_n{{0.459105,0,0,0.1,0.885347,0,0,0.1}};
    const std::array<double,8> obtained_projected_values_n{{v_n_29[0], v_n_29[1], v_n_29[2], p_n_29, v_n_53[0], v_n_53[1], v_n_53[2], p_n_53}};
    for (std::size_t i = 0; i < 8; ++i) {
        KRATOS_EXPECT_NEAR(obtained_projected_values_n[i], expected_projected_values_n[i], tol);
    }
}

}  // namespace Kratos::Testing.
