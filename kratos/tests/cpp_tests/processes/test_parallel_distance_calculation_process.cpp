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

// Project includes
#include "containers/model.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/expect.h"
// #include "includes/gid_io.h" // Include this for debugging
#include "processes/structured_mesh_generator_process.h"
#include "processes/parallel_distance_calculation_process.h"
#include "testing/testing.h"

namespace Kratos ::Testing {

namespace ParallelDistanceCalculationProcessTestInternals
{
    void SetUpDistanceField(
        ModelPart& rModelPart,
        std::function<double(Node& rNode)>& rDistanceFunction,
        std::function<double&(Node& rNode)>& rDistanceGetter)
    {
        // Set the intersected elements
        // First set an auxiliary level set field
        for (auto& r_node : rModelPart.Nodes()) {
            rDistanceGetter(r_node) = rDistanceFunction(r_node);
        }

        // Tag the nodes belonging to an intersected element0
        std::size_t n_pos, n_neg;
        for (auto& r_element : rModelPart.Elements()) {
            // Check if the element is split
            n_pos = 0;
            n_neg = 0;
            for (auto& r_node : r_element.GetGeometry()) {
                if (rDistanceGetter(r_node) > 0) {
                    n_pos++;
                } else {
                    n_neg++;
                }
            }
            // If the element is split tag its nodes
            if (n_pos != 0 && n_neg != 0) {
                for (auto& r_node : r_element.GetGeometry()) {
                    r_node.Set(SELECTED, true);
                }
            }
        }

        // Remove the nodal values from the non-intersected elements nodes and reset the flag
        for (auto& r_node : rModelPart.Nodes()) {
            if (!r_node.Is(SELECTED)) {
                double& r_dist = rDistanceGetter(r_node);
                r_dist = r_dist < 0.0 ? -1.0 : 1.0;
            } else {
                r_node.Set(SELECTED, false);
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(ParallelDistanceProcessQuadrilateral2D, KratosCoreFastSuite)
{
    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 10.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 10.0, 10.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 10.0, 0.0, 0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions" : 7,
        "element_name" : "Element2D3N",
        "create_skin_sub_model_part" : false
    })");
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("MainModelPart");
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    // Set up the intersected elements distance
    std::function<double(Node& rNode)> nodal_value_function = [](Node& rNode){return rNode.X() + rNode.Y() - 100.0/9.9;};
    std::function<double&(Node& rNode)> distance_getter = [](Node& rNode)->double&{return rNode.FastGetSolutionStepValue(DISTANCE);};
    ParallelDistanceCalculationProcessTestInternals::SetUpDistanceField(r_model_part, nodal_value_function, distance_getter);

    // Compute distance
    Parameters parallel_distance_settings(R"({
        "model_part_name" : "MainModelPart",
        "distance_variable" : "DISTANCE",
        "nodal_area_variable" : "NODAL_AREA",
        "distance_database" : "nodal_historical",
        "max_levels" : 10.0,
        "max_distance" : 10.0,
        "calculate_exact_distances_to_plane" : false
    })");
    ParallelDistanceCalculationProcess<2>(r_model_part, parallel_distance_settings).Execute();

    // // Auxiliary output for debugging
    // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/main_model_part", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_fluid.InitializeMesh(0.00);
    // gid_io_fluid.WriteMesh(r_model_part.GetMesh());
    // gid_io_fluid.FinalizeMesh();
    // gid_io_fluid.InitializeResults(0, r_model_part.GetMesh());
    // gid_io_fluid.WriteNodalResults(DISTANCE, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.WriteNodalResults(NODAL_AREA, r_model_part.Nodes(), 0, 0);
    // gid_io_fluid.FinalizeResults();

    // Check results
    const double tolerance = 1.0e-8;
    const std::array<std::size_t,4> nodal_ids = {1,28,37,64};
    const std::array<double, 4> exact_dist = {-7.14249273926, -1.08157747194, 0.93872761716, 6.99964288447};
    for (std::size_t i = 0; i < nodal_ids.size(); ++i) {
        const auto& r_node = r_model_part.GetNode(nodal_ids[i]);
        const double dist = r_node.FastGetSolutionStepValue(DISTANCE);
        // std::cout << std::setprecision(12) << dist << std::endl; // Output to update test values
        KRATOS_EXPECT_NEAR(dist, exact_dist[i], tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ParallelDistanceProcessQuadrilateralNonHistorical2D, KratosCoreFastSuite)
{
    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 10.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 10.0, 10.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 10.0, 0.0, 0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions" : 7,
        "element_name" : "Element2D3N",
        "create_skin_sub_model_part" : false
    })");
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("MainModelPart");
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    // Set up the intersected elements distance
    std::function<double(Node& rNode)> nodal_value_function = [](Node& rNode){return rNode.X() + rNode.Y() - 100.0/9.9;};
    std::function<double&(Node& rNode)> distance_getter = [](Node& rNode)->double&{return rNode.GetValue(DISTANCE);};
    ParallelDistanceCalculationProcessTestInternals::SetUpDistanceField(r_model_part, nodal_value_function, distance_getter);

    // Compute distance
    Parameters parallel_distance_settings(R"({
        "model_part_name" : "MainModelPart",
        "distance_variable" : "DISTANCE",
        "nodal_area_variable" : "NODAL_AREA",
        "distance_database" : "nodal_non_historical",
        "max_levels" : 10.0,
        "max_distance" : 10.0,
        "calculate_exact_distances_to_plane" : false
    })");
    ParallelDistanceCalculationProcess<2>(r_model_part, parallel_distance_settings).Execute();

    // Check results
    const double tolerance = 1.0e-8;
    const std::array<std::size_t,4> nodal_ids = {1,28,37,64};
    const std::array<double, 4> exact_dist = {-7.14249273926, -1.08157747194, 0.93872761716, 6.99964288447};
    for (std::size_t i = 0; i < nodal_ids.size(); ++i) {
        const auto& r_node = r_model_part.GetNode(nodal_ids[i]);
        const double dist = r_node.GetValue(DISTANCE);
        KRATOS_EXPECT_NEAR(dist, exact_dist[i], tolerance);
    }
}

}  // namespace Kratos::Testing.
