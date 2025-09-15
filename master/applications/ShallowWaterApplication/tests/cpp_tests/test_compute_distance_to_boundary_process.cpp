//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/variable_utils.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "custom_processes/calculate_distance_to_boundary_process.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::NodeType NodeType;

void InitializeAndFillModelPart(ModelPart& rModelPart)
{
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  5,
        "element_name":         "Element2D3N",
        "condition_name":       "LineCondition"
    })");

    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDistanceToOneBoundaryProcess, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("model_part");
    InitializeAndFillModelPart(r_model_part);

    auto& r_boundary_part = r_model_part.CreateSubModelPart("boundary");
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.X() >= 1.0) {
            r_boundary_part.AddNode(&r_node);
        }
    }

    CalculateDistanceToBoundaryProcess(r_model_part, r_boundary_part).ExecuteBeforeSolutionLoop();

    const double tolerance = 1e-16;
    for (auto& r_node : r_model_part.Nodes()) {
        const double distance = 1.0 - r_node.X();
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(DISTANCE), distance, tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDistanceToTwoBoundariesProcess, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("model_part");
    InitializeAndFillModelPart(r_model_part);

    auto& r_boundary_part_1 = r_model_part.CreateSubModelPart("boundary_1");
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.X() >= 1.0) {
            r_boundary_part_1.AddNode(&r_node);
        }
    }

    auto& r_boundary_part_2 = r_model_part.CreateSubModelPart("boundary_2");
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.Y() >= 1.0) {
            r_boundary_part_2.AddNode(&r_node);
        }
    }

    CalculateDistanceToBoundaryProcess proc_1(r_model_part, r_boundary_part_1);
    CalculateDistanceToBoundaryProcess proc_2(r_model_part, r_boundary_part_2);
    proc_1.ExecuteBeforeSolutionLoop();
    proc_2.ExecuteBeforeSolutionLoop();

    const double tolerance = 1e-16;
    for (auto& r_node : r_model_part.Nodes()) {
        const double distance = std::min(1.0 - r_node.X(), 1.0 - r_node.Y());
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(DISTANCE), distance, tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDistanceToTwoBoundariesBruteForceProcess, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("model_part");
    InitializeAndFillModelPart(r_model_part);

    auto& r_boundary_part = r_model_part.CreateSubModelPart("boundary");
    for (auto& r_node : r_model_part.Nodes()) {
        if (r_node.X() >= 1.0) {
            r_boundary_part.AddNode(&r_node);
        }
        if ((r_node.Y() >= 1.0)) {
            r_boundary_part.AddNode(&r_node);
        }
    }

    CalculateDistanceToBoundaryProcess(r_model_part, r_boundary_part).ExecuteBeforeSolutionLoop();

    const double tolerance = 1e-16;
    for (auto& r_node : r_model_part.Nodes()) {
        const double distance = std::min(1.0 - r_node.X(), 1.0 - r_node.Y());
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(DISTANCE), distance, tolerance);
    }
}

} // namespace Testing

} // namespace Kratos
