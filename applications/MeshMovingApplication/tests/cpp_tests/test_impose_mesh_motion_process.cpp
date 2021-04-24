//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// Project includes
#include "containers/model.h"
#include "custom_processes/impose_mesh_motion_process.h"
#include "testing/testing.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/mesh_moving_variables.h"

namespace Kratos
{
namespace Testing
{


void FillModelPart(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);

    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

    // Generate mesh
    using GeometryType = Quadrilateral2D4<Node<3>>;

    GeometryType::PointsArrayType domain_corners;
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(4, 0.0, 2.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(3, 2.0, 2.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(2, 2.0, 0.0, 0.0));

    Parameters mesh_generator_parameters(R"(
    {
        "number_of_divisions" : 2,
        "element_name" : "Element2D3N",
        "create_skin_sub_model_part" : false
    })");

    StructuredMeshGeneratorProcess(
        GeometryType(domain_corners),
        rModelPart,
        mesh_generator_parameters
    ).Execute();

    // Add DoFs
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }
}


KRATOS_TEST_CASE_IN_SUITE(ImposeMeshMotionProcessRotationAxis, MeshMovingApplicationFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 3);
    FillModelPart(r_model_part);

    // Rotate around an offset y-axis, then translate by [1,2,3]
    Parameters parameters(R"(
    {
        "rotation_definition" : "rotation_axis",
        "rotation_axis"       : [0, 1, 0],
        "reference_point"     : [-1, 0, 0],
        "rotation_angle"      : -1.57079632679,
        "translation_vector"  : [1, 2, 3]
    })");

    ImposeMeshMotionProcess impose_process(r_model_part, parameters);
    impose_process.ExecuteInitializeSolutionStep();

    Node<3> node;
    const double tolerance = 1e-10;

    node = r_model_part.GetNode(1);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(2);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(3);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(4);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(5);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(6);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(7);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);

    node = r_model_part.GetNode(8);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);

    node = r_model_part.GetNode(9);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);
}


KRATOS_TEST_CASE_IN_SUITE(ImposeMeshMotionProcessEulerAngles, MeshMovingApplicationFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 3);
    FillModelPart(r_model_part);

    // Rotate around an offset y-axis, then translate by [1,2,3]
    Parameters parameters(R"(
    {
        "rotation_definition" : "euler_angles",
        "euler_angles"        : [-1.57079632679, -1.57079632679, 1.57079632679],
        "reference_point"     : [-1, 0, 0],
        "translation_vector"  : [1, 2, 3]
    })");

    ImposeMeshMotionProcess impose_process(r_model_part, parameters);
    impose_process.ExecuteInitializeSolutionStep();

    Node<3> node;
    const double tolerance = 1e-10;

    node = r_model_part.GetNode(1);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(2);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(3);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), 0.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 4.00, tolerance);

    node = r_model_part.GetNode(4);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(5);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(6);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -1.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 5.00, tolerance);

    node = r_model_part.GetNode(7);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);

    node = r_model_part.GetNode(8);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);

    node = r_model_part.GetNode(9);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_X), -2.0, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Y), 2.00, tolerance);
    KRATOS_CHECK_NEAR(node.GetValue(MESH_DISPLACEMENT_Z), 6.00, tolerance);
}


} // namespace Testing
} // namespace Kratos