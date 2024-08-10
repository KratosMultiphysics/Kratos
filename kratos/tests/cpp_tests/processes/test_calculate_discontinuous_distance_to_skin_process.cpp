//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//

// Project includes
#include "testing/testing.h"
#include "intrusive_ptr/intrusive_ptr.hpp"
#include "containers/model.h"
#include "includes/expect.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessHorizontalPlane2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.AddNodalSolutionStepVariable(DISTANCE);
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.5;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
    skin_part.CreateNewNode(2,  1.0, plane_height, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.5,-0.5,0.5};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_edge_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.5,0.5,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximation2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.AddNodalSolutionStepVariable(DISTANCE);
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.5;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
    skin_part.CreateNewNode(2, 0.75, plane_height, 0.0);
    skin_part.CreateNewNode(3, 0.75, -1.0, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 2, {{2, 3}}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.483157,0.0216888,0.380052};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_edge_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.375,0.5,0.75};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessCubeInCube3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, -0.5);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2,  0.5, -0.5, -0.5);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, -0.5);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, -0.5,  0.5, -0.5);
    Node::Pointer p_point_5 = Kratos::make_intrusive<Node>(5, -0.5, -0.5,  0.5);
    Node::Pointer p_point_6 = Kratos::make_intrusive<Node>(6,  0.5, -0.5,  0.5);
    Node::Pointer p_point_7 = Kratos::make_intrusive<Node>(7,  0.5,  0.5,  0.5);
    Node::Pointer p_point_8 = Kratos::make_intrusive<Node>(8, -0.5,  0.5,  0.5);

    Hexahedra3D8<Node > geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   5,
        "element_name":     "Element3D4N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(VELOCITY);
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the cube skin
    const double cube_radious = 0.35;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(1, -cube_radious, -cube_radious, -cube_radious);
    skin_part.CreateNewNode(2,  cube_radious, -cube_radious, -cube_radious);
    skin_part.CreateNewNode(3,  cube_radious,  cube_radious, -cube_radious);
    skin_part.CreateNewNode(4, -cube_radious,  cube_radious, -cube_radious);
    skin_part.CreateNewNode(5, -cube_radious, -cube_radious,  cube_radious);
    skin_part.CreateNewNode(6,  cube_radious, -cube_radious,  cube_radious);
    skin_part.CreateNewNode(7,  cube_radious,  cube_radious,  cube_radious);
    skin_part.CreateNewNode(8, -cube_radious,  cube_radious,  cube_radious);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
    skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    parameters.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_dist_elem_1 = (volume_part.ElementsBegin() + 7)->GetValue(ELEMENTAL_DISTANCES);
    const auto &r_dist_elem_2 = (volume_part.ElementsEnd() - 7)->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values_elem_1 = {-0.15,0.05,0.05,-0.15};
    const std::vector<double> expected_values_elem_2 = {-0.05,0.15,0.15,-0.05};
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_1, expected_values_elem_1, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_2, expected_values_elem_2, 1.0e-6);

    // Check edge distances
    const auto &r_dist_elem_1_edge = (volume_part.ElementsBegin() + 7)->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_dist_elem_2_edge = (volume_part.ElementsEnd() - 7)->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_elem_1_edge = {0.75,-1.0,0.25,-1.0,0.25,0.25};
    const std::vector<double> expected_values_elem_2_edge = {0.25,-1.0,0.75,-1.0,0.75,0.75};
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_1_edge, expected_values_elem_1_edge, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_2_edge, expected_values_elem_2_edge, 1.0e-6);

    //Check extra edge distances - elem_1 has 4 cut edges; elem_2 has only 3 cut edges, but is still completely intersected
    const auto &r_dist_elem_1_edge_extra = (volume_part.ElementsBegin() + 7)->GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_dist_elem_2_edge_extra = (volume_part.ElementsEnd() - 2)->GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_values_elem_1_edge_extra = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_2_edge_extra = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_1_edge_extra, expected_values_elem_1_edge_extra, 1.0e-10);
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_elem_2_edge_extra, expected_values_elem_2_edge_extra, 1.0e-10);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSharpCornerInCube3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0);
    Node::Pointer p_point_5 = Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0);
    Node::Pointer p_point_6 = Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0);
    Node::Pointer p_point_7 = Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0);
    Node::Pointer p_point_8 = Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0);

    Hexahedra3D8<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   1,
        "element_name":     "Element3D4N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the corner entities
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(1, 2.0, 0.5, 2.0);
    skin_part.CreateNewNode(2, 2.0, 0.5, -2.0);
    skin_part.CreateNewNode(3, 0.25, 0.5, 2.0);
    skin_part.CreateNewNode(4, 0.25, 0.5, -2.0);
    skin_part.CreateNewNode(5, 0.25, -2.0, 2.0);
    skin_part.CreateNewNode(6, 0.25, -2.0, -2.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties);
    skin_part.CreateNewElement("Element3D3N", 2, {2, 3, 4}, p_properties);
    skin_part.CreateNewElement("Element3D3N", 3, {3, 4, 5}, p_properties);
    skin_part.CreateNewElement("Element3D3N", 4, {4, 5, 6}, p_properties);

    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_dist_begin = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const auto &r_dist_end = (volume_part.ElementsEnd() - 1)->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values_begin = {1.73205,1.73205,1.73205,1.73205};
    const std::vector<double> expected_values_end = {-0.406059,-0.489839,0.388306,0.0649427};
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_begin, expected_values_begin, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_end, expected_values_end, 1.0e-6);

    // Check edge distances
    const auto &r_dist_begin_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_dist_end_edge = (volume_part.ElementsEnd() - 1)->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_begin_edge = {-1,-1,0.625,0.375,-1,-1};
    const std::vector<double> expected_values_end_edge = {-1,0.5,0.5,0.75,0.75,0.625};
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_begin_edge, expected_values_begin_edge, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_dist_end_edge, expected_values_end_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSinglePointTangent2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    const double plane_height = 1.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
    skin_part.CreateNewNode(2,  1.0, plane_height, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    KRATOS_EXPECT_NEAR(r_elem_dist[0], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSinglePointTangentOnNodeInteresected2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, 1.0, 0.0);
    skin_part.CreateNewNode(2,  0.0, 1.0, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const auto &r_to_split = (fluid_part.ElementsBegin())->Is(TO_SPLIT);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_TRUE(r_to_split);
    KRATOS_EXPECT_NEAR(r_elem_dist[0], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSinglePointTangentOnNodeNotInteresected2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, 1.0, 0.0);
    skin_part.CreateNewNode(2,  0.0, 1.0, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    Parameters parameters;
    parameters.AddBool("use_positive_epsilon_for_zero_values", false);
    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const auto &r_to_split = (fluid_part.ElementsBegin())->Is(TO_SPLIT);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_TRUE(!r_to_split);
    KRATOS_EXPECT_NEAR(r_elem_dist[0], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], -epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSingleLineTangent2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
    skin_part.CreateNewNode(2,  1.0, plane_height, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_NEAR(r_elem_dist[0], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], 1.0, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSinglePointAndManyIntersectEdge2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, 1.0, 0.0);
    skin_part.CreateNewNode(2,  1.0, 1.0, 0.0);
    skin_part.CreateNewNode(3,  1.0/6.0, 2.0/3.0, 0.0);
    skin_part.CreateNewNode(4,  1.0, 2.0/3.0, 0.0);
    skin_part.CreateNewNode(5,  2.0/3.0, 1.0/6.0, 0.0);
    skin_part.CreateNewNode(6,  1.0, 1.0/6.0, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 2, {{2, 3}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 3, {{3, 4}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 4, {{4, 5}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 5, {{5, 6}}, p_properties_1);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_NEAR(r_elem_dist[0], -std::sqrt(2)/2.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessMultipleTangent2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 2.0, -1.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 1.0, 0.0);
    fluid_part.CreateNewNode(3, 1.2, -1.0, 0.0);
    fluid_part.CreateNewNode(4, 2.6, 0.0, 0.0);
    fluid_part.CreateNewNode(5, 0.3, -0.4, 0.0);
    fluid_part.CreateNewNode(6, 1.5, 0.0, 0.0);
    fluid_part.CreateNewNode(7, 2.0, 1.0, 0.0);
    fluid_part.CreateNewNode(8, 0.5, 0.5, 0.0);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element2D3N", 1, {1, 4, 6}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 2, {3, 6, 5}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 3, {6, 4, 7}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 4, {1, 6, 3}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 5, {6, 2, 8}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 6, {6, 7, 2}, p_properties_0);
    fluid_part.CreateNewElement("Element2D3N", 7, {6, 8, 5}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -3.0, plane_height, 0.0);
    skin_part.CreateNewNode(2,  3.0, plane_height, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    for (auto& r_node : fluid_part.Nodes()) {
        r_node.SetValue(DISTANCE, std::numeric_limits<double>::max());
    }
    for (auto& r_elem : fluid_part.Elements()) {
        const auto& r_elem_dist = r_elem.GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < r_elem_dist.size(); i++) {
            auto& r_nodal_dist = r_elem.GetGeometry()[i].GetValue(DISTANCE);
            if (r_elem_dist[i] < r_nodal_dist) {
                r_elem.GetGeometry()[i].SetValue(DISTANCE, r_elem_dist[i]);
            }
        }
    }

    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    // Check values
    for (auto& r_node : fluid_part.Nodes()) {
        if (std::abs(r_node.Y()) < epsilon) {
            KRATOS_EXPECT_NEAR(r_node.GetValue(DISTANCE), epsilon, 1e-16);
        } else {
            KRATOS_EXPECT_NEAR(r_node.GetValue(DISTANCE), r_node.Y(), 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSinglePointTangent3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0 , 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    fluid_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin
    const double plane_height = 1.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, -5.0, 5.0, plane_height);
    skin_part.CreateNewNode(902, -5.0, -5.0, plane_height);
    skin_part.CreateNewNode(903, 5.0, -5.0, plane_height);
    skin_part.CreateNewNode(904, 5.0, 5.0, plane_height);
    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);


    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_NEAR(r_elem_dist[0], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], -1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[3], epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSingleLineTangent3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0 , 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    fluid_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, 0.0, -5.0, 5.0);
    skin_part.CreateNewNode(902, 0.0, -5.0, -5.0);
    skin_part.CreateNewNode(903, 0.0, 5.0, -5.0);
    skin_part.CreateNewNode(904, 0.0, 5.0, 5.0);
    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);


    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_NEAR(r_elem_dist[0], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], 1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], 1.0, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[3], epsilon, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSingleFaceTangent3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0 , 0.0, 0.0);
    fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    fluid_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    fluid_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, -5.0, 5.0, plane_height);
    skin_part.CreateNewNode(902, -5.0, -5.0, plane_height);
    skin_part.CreateNewNode(903, 5.0, -5.0, plane_height);
    skin_part.CreateNewNode(904, 5.0, 5.0, plane_height);
    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);


    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    KRATOS_EXPECT_NEAR(r_elem_dist[0], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[1], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[2], epsilon, 1e-16);
    KRATOS_EXPECT_NEAR(r_elem_dist[3], 1.0, 1e-16);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessMultipleTangent3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.0 , 0.0, 0.0);
    fluid_part.CreateNewNode(2, 0.0, 1.0, 1.0);
    fluid_part.CreateNewNode(3, 1.0, -1.0, 1.0);
    fluid_part.CreateNewNode(4, -2.0, 1.0, 1.0);
    fluid_part.CreateNewNode(5, 0.0, 1.0, -1.0);
    fluid_part.CreateNewNode(6, 1.0, -1.0, -1.0);
    fluid_part.CreateNewNode(7, -2.0, 1.0, -1.0);
    fluid_part.CreateNewNode(8, 1.0, 1.0, 0.0);
    fluid_part.CreateNewNode(9, -2.0, -0.5, 0.0);
    fluid_part.CreateNewNode(10, -3.0, 1.0, 0.0);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element3D4N", 1, {1, 3, 2, 4}, p_properties_0);
    fluid_part.CreateNewElement("Element3D4N", 2, {1, 5, 6, 7}, p_properties_0);
    fluid_part.CreateNewElement("Element3D4N", 3, {1, 6, 5, 8}, p_properties_0);
    fluid_part.CreateNewElement("Element3D4N", 4, {1, 9, 10, 7}, p_properties_0);

    // Generate the skin
    const double plane_height = 0.0;
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, -5.0, 5.0, plane_height);
    skin_part.CreateNewNode(902, -5.0, -5.0, plane_height);
    skin_part.CreateNewNode(903, 5.0, -5.0, plane_height);
    skin_part.CreateNewNode(904, 5.0, 5.0, plane_height);
    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();

    // Check values
    const double epsilon = std::numeric_limits<double>::epsilon()*1e3;
    for (auto &elem : fluid_part.Elements()) {
        const auto& r_elem_dist = elem.GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < r_elem_dist.size(); i++) {
            if (std::abs(elem.GetGeometry()[i].Z()) < epsilon) {
                KRATOS_EXPECT_NEAR(r_elem_dist[i], epsilon, 1e-16);
            } else {
                KRATOS_EXPECT_NEAR(r_elem_dist[i], elem.GetGeometry()[i].Z(), 1e-16);
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessHorizontalPlane3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the evil element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.214286, -0.357143, 0.0714286);
    volume_part.CreateNewNode(2, 0.214286, -0.214286, 0.0714286);
    volume_part.CreateNewNode(3, 0.357143, -0.214286, 0.0714286);
    volume_part.CreateNewNode(4, 0.214286, -0.357143, -0.0714286);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the cube skin
    const double plane_height = 0.0;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, -1.0, plane_height);
    skin_part.CreateNewNode(2,  1.0, -1.0, plane_height);
    skin_part.CreateNewNode(3,  1.0,  1.0, plane_height);
    skin_part.CreateNewNode(4, -1.0,  1.0, plane_height);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {1, 2, 4}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {2, 3, 4}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {0.0714286,0.0714286,0.0714286,-0.0714286};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {-1,-1,-1,0.5,0.5,0.5};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-6);
}



KRATOS_TEST_CASE_IN_SUITE(HorizontalPlaneZeroDiscontinuousDistanceProcess, KratosCoreFastSuite)
{

    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point1 = Kratos::make_intrusive<Node>(1, 0.00, 0.00, 0.00);
    Node::Pointer p_point2 = Kratos::make_intrusive<Node>(2, 10.00, 0.00, 0.00);
    Node::Pointer p_point3 = Kratos::make_intrusive<Node>(3, 10.00, 10.00, 0.00);
    Node::Pointer p_point4 = Kratos::make_intrusive<Node>(4, 0.00, 10.00, 0.00);
    Node::Pointer p_point5 = Kratos::make_intrusive<Node>(5, 0.00, 0.00, 10.00);
    Node::Pointer p_point6 = Kratos::make_intrusive<Node>(6, 10.00, 0.00, 10.00);
    Node::Pointer p_point7 = Kratos::make_intrusive<Node>(7, 10.00, 10.00, 10.00);
    Node::Pointer p_point8 = Kratos::make_intrusive<Node>(8, 0.00, 10.00, 10.00);

    Hexahedra3D8<Node > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   2,
        "element_name":     "Element3D4N"
    })");

    Model current_model;
    ModelPart &volume_part = current_model.CreateModelPart("Volume");
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, 0.0, 0.0, 5.0);
    skin_part.CreateNewNode(902, 10.0, 0.0, 5.0);
    skin_part.CreateNewNode(903, 10.0, 10.0, 5.0);
    skin_part.CreateNewNode(904, 0.0, 10.0, 5.0);
    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);

    // Compute distance
    CalculateDiscontinuousDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

    // Check values
    for (auto &elem : volume_part.Elements()) {
        const auto& r_elem_dist = elem.GetValue(ELEMENTAL_DISTANCES);
        for (unsigned int i = 0; i < r_elem_dist.size(); i++) {
            double distance = elem.GetGeometry()[i].Z()-5.0;
            KRATOS_EXPECT_NEAR(r_elem_dist[i], distance, 1e-6);
        }
    }

}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationSkewed3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the tetrahedron element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.0, -0.5, 0.0);
    volume_part.CreateNewNode(2, 1.0, -0.5, 0.0);
    volume_part.CreateNewNode(3, 0.0,  0.5, 0.0);
    volume_part.CreateNewNode(4, 0.0, -0.5, 1.0);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin such that one edge is cut twice
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -1.0, -1.0,  0.75);
    skin_part.CreateNewNode(2,  1.0, -1.0,  0.75);
    skin_part.CreateNewNode(3, -1.0,  1.0,  0.75);
    skin_part.CreateNewNode(4, 0.75, -1.0,  1.0);
    skin_part.CreateNewNode(5, 0.75, -1.0, -1.0);
    skin_part.CreateNewNode(6, 0.75,  1.0, -1.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {1,2,3}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {4,5,6}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.569400,0.1044875,-0.266495,0.104487};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.75,0.25,-1,0.75,0.5,0.75};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationVertical3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the triangular element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    volume_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin such that there is 4 intersection pts.
    // Recall that with more than 3 intersection pts. the plane
    // approximation is used. Since the skin in here yields a
    // uniplanar intersection, the approximated plane is the
    // same one as the original intersection one.
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(1, 0.5, -1.0,  1.0);
    skin_part.CreateNewNode(2, 0.5, -1.0, -1.0);
    skin_part.CreateNewNode(3, 0.5,  1.0, -1.0);
    skin_part.CreateNewNode(4, 0.5,  1.0,  1.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {1,2,4}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {2,3,4}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.5,0.5,-0.5,-0.5};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-10);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.5,0.5,-1,-1,0.5,-1};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessOneEdgeIntersection3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the tetrahedron element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.666963, 0.800762, 0.388769);
    volume_part.CreateNewNode(2, 0.731067, 0.821936, 0.422077);
    volume_part.CreateNewNode(3, 0.652002, 0.85453, 0.463652);
    volume_part.CreateNewNode(4, 0.684484, 0.796908, 0.48275);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin such that it only intersects in one edge
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, 0.675, 0.803109, 0.5);
    skin_part.CreateNewNode(2, 0.663088, 0.808771, 0.476277);
    skin_part.CreateNewNode(3, 0.685008, 0.796367, 0.479053);
    skin_part.CreateNewNode(4, 0.682845, 0.794215, 0.449949);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {1,2,3}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {3,2,4}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {0.135661,0.135661,0.135661,0.135661};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-5);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {-1,-1,-1,0.959824,-1,-1};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessMultipleIntersections3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the tetrahedron element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.597905, 0.597905, 0.100929);
    volume_part.CreateNewNode(2, 0.608229, 0.490745, 0.204129);
    volume_part.CreateNewNode(3, 0.697865, 0.493815, 0.126859);
    volume_part.CreateNewNode(4, 0.65897, 0.571858, 0.175526);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin such that it has multiple intersections
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, 0.633131, 0.539808, 0.178766);
    skin_part.CreateNewNode(2, 0.671961, 0.517362, 0.195651);
    skin_part.CreateNewNode(3, 0.66866, 0.566563, 0.200629);
    skin_part.CreateNewNode(4, 0.635672, 0.588229, 0.189664);
    skin_part.CreateNewNode(5, 0.631307, 0.501763, 0.175569);
    skin_part.CreateNewNode(6, 0.664311, 0.467496, 0.19268);
    skin_part.CreateNewNode(7, 0.595066, 0.567432, 0.169977);
    skin_part.CreateNewNode(8, 0.589401, 0.523455, 0.162424);
    skin_part.CreateNewNode(9, 0.629759, 0.452428, 0.178442);
    skin_part.CreateNewNode(10, 0.591444, 0.479469, 0.162781);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {1,2,3}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {1,3,4}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 3, {1,5,2}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 4, {2,5,6}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 5, {7,1,4}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 6, {7,8,1}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 7, {8,5,1}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 8, {9,6,5}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 9, {9,10,5}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 10, {5,10,9}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distance values
    const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.0636738,0.0342287,-0.0709816,-0.0159295};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.65038,0.325336,-1,-1,0.682415,-1};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessStandard3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the tetrahedron element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, -0.625, -1.625, 1.125);
    volume_part.CreateNewNode(2, -0.5, -1.75, 1.25);
    volume_part.CreateNewNode(3, -0.5, -1.75, 1.0);
    volume_part.CreateNewNode(4, -0.375, -1.625, 1.125);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin such that it generates a standard intersection
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(2196, -0.542249, -1.7162, 1.23726);
    skin_part.CreateNewNode(2155, -0.544766, -1.71316, 1.19334);
    skin_part.CreateNewNode(2228, -0.507166, -1.69137, 1.20104);
    skin_part.CreateNewNode(2220, -0.502438, -1.68588, 1.17549);
    skin_part.CreateNewNode(2335, -0.457592, -1.66163, 1.20154);
    skin_part.CreateNewNode(2118, -0.547284, -1.71012, 1.14943);
    skin_part.CreateNewNode(2368, -0.419474, -1.62945, 1.14036);
    skin_part.CreateNewNode(2318, -0.453428, -1.65552, 1.17619);
    skin_part.CreateNewNode(2276, -0.462077, -1.65634, 1.14339);
    skin_part.CreateNewNode(2324, -0.42378, -1.62452, 1.09704);
    skin_part.CreateNewNode(2422, -0.389158, -1.5991, 1.10588);
    skin_part.CreateNewNode(2263, -0.463227, -1.65498, 1.12885);
    skin_part.CreateNewNode(2250, -0.464376, -1.65362, 1.11431);
    skin_part.CreateNewNode(2236, -0.460832, -1.64688, 1.08915);
    skin_part.CreateNewNode(2209, -0.503559, -1.68456, 1.16095);
    skin_part.CreateNewNode(2190, -0.50468, -1.68323, 1.14641);
    skin_part.CreateNewNode(2086, -0.550732, -1.70604, 1.10582);
    skin_part.CreateNewNode(2154, -0.513394, -1.68396, 1.11361);
    skin_part.CreateNewNode(2141, -0.509278, -1.67779, 1.08826);
    skin_part.CreateNewNode(2049, -0.55418, -1.70196, 1.06221);
    skin_part.CreateNewNode(2131, -0.510714, -1.67615, 1.07382);
    skin_part.CreateNewNode(2120, -0.512149, -1.6745, 1.05938);
    skin_part.CreateNewNode(2091, -0.521568, -1.67444, 1.02683);
    skin_part.CreateNewNode(2021, -0.558662, -1.69672, 1.01895);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 32933, {2196,2155,2228}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 32940, {2228,2220,2335}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 32934, {2228,2118,2220}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 32937, {2155,2118,2228}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38354, {2368,2318,2276}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33145, {2368,2324,2422}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33157, {2368,2276,2263}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33163, {2368,2263,2250}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38573, {2236,2368,2250}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38576, {2324,2368,2236}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38353, {2209,2318,2220}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38355, {2318,2335,2220}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38357, {2209,2276,2318}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38356, {2118,2209,2220}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38351, {2209,2190,2276}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33155, {2118,2086,2154}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33161, {2118,2154,2190}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38350, {2118,2190,2209}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38579, {2236,2250,2141}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33160, {2263,2154,2250}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33162, {2154,2141,2250}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33164, {2263,2190,2154}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33159, {2086,2049,2154}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38577, {2131,2236,2141}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33156, {2154,2049,2141}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38580, {2049,2131,2141}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 38574, {2049,2120,2131}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33399, {2049,2091,2120}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 33393, {2049,2021,2091}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.108523,0.0485713,0.0764035,0.0222594};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {0.690813,-1,0.413157,0.829797,-1,-1};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessBoundaryIntersection3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the tetrahedron element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(787, 0.3, 0.1, 0.5);
    volume_part.CreateNewNode(629, 0.3, 0.2, 0.5);
    volume_part.CreateNewNode(700, 0.4, 0.2, 0.5);
    volume_part.CreateNewNode(712, 0.3, 0.1, 0.4);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 861, {787, 629, 700, 712}, p_properties_0);

    // Generate the skin such that one edge intersection is close to the boundary entity
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(345, 0.372131, 0.174194, 0.5);
    skin_part.CreateNewNode(375, 0.396836, 0.16555, 0.5);
    skin_part.CreateNewNode(333, 0.384461, 0.170563, 0.475061);
    skin_part.CreateNewNode(384, 0.384461, 0.170563, 0.524939);
    skin_part.CreateNewNode(351, 0.35857, 0.180817, 0.524898);
    skin_part.CreateNewNode(310, 0.348141, 0.184661, 0.5);
    skin_part.CreateNewNode(319, 0.335115, 0.192276, 0.524874);
    skin_part.CreateNewNode(280, 0.325, 0.196891, 0.5);
    skin_part.CreateNewNode(266, 0.335115, 0.192276, 0.475126);
    skin_part.CreateNewNode(290, 0.31342, 0.204924, 0.524859);
    skin_part.CreateNewNode(248, 0.302838, 0.210816, 0.5);
    skin_part.CreateNewNode(236, 0.31342, 0.204924, 0.475141);
    skin_part.CreateNewNode(299, 0.35857, 0.180817, 0.475102);
    skin_part.CreateNewNode(285, 0.371864, 0.179378, 0.442703);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 554, {345,375,333}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1081, {375,345,384}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1295, {384,345,351}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1199, {351,310,319}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 574, {280,310,266}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1101, {310,280,319}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1193, {319,280,290}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 582, {248,280,236}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1109, {280,248,290}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 667, {280,266,236}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 562, {310,345,299}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 1089, {345,310,351}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 769, {345,333,299}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 673, {310,299,266}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 770, {299,333,285}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Check elemental distances
    const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values = {-0.0984855,-0.00883326,0.0352186,-0.103167};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist, expected_values, 1.0e-6);

    // Check edge distances
    const auto &r_elem_dist_edge = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_edge = {-1,0.200519,0.263407,-1,-1,0.254497};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_edge, expected_values_edge, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessDoubleEmbeddedVariableComplex, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a triangle element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(TEMPERATURE);
    skin_part.CreateNewNode(1,-0.1, 0.5,0.0);
    skin_part.CreateNewNode(2, 0.1, 0.5,0.0);
    skin_part.CreateNewNode(3,-0.1, 0.3,0.0);
    skin_part.CreateNewNode(4, 0.1, 0.3,0.0);
    skin_part.CreateNewNode(5, 0.1,-0.1,0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 2, {{2,3}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 3, {{3,4}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 4, {{4,5}}, p_properties_1);

    // Set the embedded cube double variable
    for (auto &i_node : skin_part.Nodes()) {
        i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.Y();
    }

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part);
    disc_dist_proc.Execute();
    disc_dist_proc.CalculateEmbeddedVariableFromSkin(TEMPERATURE, TEMPERATURE);

    // Check values
    KRATOS_EXPECT_NEAR(volume_part.GetElement(1).GetValue(TEMPERATURE), 0.2, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessArrayEmbeddedVariableComplex, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a triangle element
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

    // Generate the skin
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(1,-0.1, 0.5,0.0);
    skin_part.CreateNewNode(2, 0.1, 0.5,0.0);
    skin_part.CreateNewNode(3,-0.1, 0.3,0.0);
    skin_part.CreateNewNode(4, 0.1, 0.3,0.0);
    skin_part.CreateNewNode(5, 0.1,-0.1,0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 2, {{2,3}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 3, {{3,4}}, p_properties_1);
    skin_part.CreateNewElement("Element2D2N", 4, {{4,5}}, p_properties_1);

    // Set the embedded cube double variable
    for (auto &i_node : skin_part.Nodes()) {
        i_node.FastGetSolutionStepValue(VELOCITY_X) = i_node.X();
        i_node.FastGetSolutionStepValue(VELOCITY_Y) = i_node.Y();
    }

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part);
    disc_dist_proc.Execute();
    disc_dist_proc.CalculateEmbeddedVariableFromSkin(VELOCITY, EMBEDDED_VELOCITY);

    // Check values
    KRATOS_EXPECT_NEAR(volume_part.GetElement(1).GetValue(EMBEDDED_VELOCITY_X), 0.05, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(1).GetValue(EMBEDDED_VELOCITY_Y), 0.2, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessDoubleEmbeddedVariable, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 7,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.AddNodalSolutionStepVariable(TEMPERATURE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the cube skin
    const double cube_radious = 0.25;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(TEMPERATURE);
    skin_part.CreateNewNode(1, -cube_radious, -cube_radious, 0.0);
    skin_part.CreateNewNode(2, -cube_radious,  cube_radious, 0.0);
    skin_part.CreateNewNode(3,  cube_radious,  cube_radious, 0.0);
    skin_part.CreateNewNode(4,  cube_radious, -cube_radious, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N",  1, {{1,2}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  2, {{2,3}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  3, {{3,4}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  4, {{4,1}}, p_properties);

    // Set the embedded cube double variable
    for (auto &i_node : skin_part.Nodes()) {
        i_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0;
    }

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part);
    disc_dist_proc.Execute();
    disc_dist_proc.CalculateEmbeddedVariableFromSkin(TEMPERATURE, TEMPERATURE);

    // Check values
    KRATOS_EXPECT_NEAR(volume_part.GetElement(16).GetValue(TEMPERATURE), 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(17).GetValue(TEMPERATURE), 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(68).GetValue(TEMPERATURE), 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(69).GetValue(TEMPERATURE), 0.0, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessArrayEmbeddedVariable, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 7,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the cube skin
    const double cube_radious = 0.25;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(1, -cube_radious, -cube_radious, 0.0);
    skin_part.CreateNewNode(2, -cube_radious,  cube_radious, 0.0);
    skin_part.CreateNewNode(3,  cube_radious,  cube_radious, 0.0);
    skin_part.CreateNewNode(4,  cube_radious, -cube_radious, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N",  1, {{1,2}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  2, {{2,3}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  3, {{3,4}}, p_properties);
    skin_part.CreateNewElement("Element2D2N",  4, {{4,1}}, p_properties);

    // Set the embedded cube array variable
    array_1d<double,3> velocity = ZeroVector(3);
    velocity[0] = 1.0;
    velocity[1] = 1.0;
    for (auto &i_node : skin_part.Nodes()) {
        i_node.FastGetSolutionStepValue(VELOCITY) = velocity;
    }

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part);
    disc_dist_proc.Execute();
    disc_dist_proc.CalculateEmbeddedVariableFromSkin(VELOCITY, EMBEDDED_VELOCITY);

    // Check values
    KRATOS_EXPECT_NEAR(volume_part.GetElement(16).GetValue(EMBEDDED_VELOCITY)[0], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(16).GetValue(EMBEDDED_VELOCITY)[1], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(16).GetValue(EMBEDDED_VELOCITY)[2], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(17).GetValue(EMBEDDED_VELOCITY)[0], 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(17).GetValue(EMBEDDED_VELOCITY)[1], 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(17).GetValue(EMBEDDED_VELOCITY)[2], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(68).GetValue(EMBEDDED_VELOCITY)[0], 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(68).GetValue(EMBEDDED_VELOCITY)[1], 1.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(68).GetValue(EMBEDDED_VELOCITY)[2], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(69).GetValue(EMBEDDED_VELOCITY)[0], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(69).GetValue(EMBEDDED_VELOCITY)[1], 0.0, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetElement(69).GetValue(EMBEDDED_VELOCITY)[2], 0.0, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessIncisedVsIntersected2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    const double plane_height = 0.25;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.4, plane_height, 0.0);
    skin_part.CreateNewNode(2,  0.4, plane_height, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    parameters.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Count intersected and incised elements
    const uint8_t n_elements = 8;
    const uint8_t n_edges = 3;
    uint8_t n_cut_edges = 0;
    uint8_t n_incised = 0;
    uint8_t n_intersected = 0;
    for (uint8_t i = 0; i < n_elements; ++i) {
        const auto &r_edge_dist = (volume_part.ElementsBegin() + i)->GetValue(ELEMENTAL_EDGE_DISTANCES);
        n_cut_edges = 0;
        for (uint8_t j = 0; j < n_edges; ++j) {
            if (r_edge_dist[j] >= 0){
                n_cut_edges++;
            }
        }
        if (n_cut_edges > 0){
            if (n_cut_edges > 1){
                n_intersected++;
            } else {
                n_incised++;
            }
        }
    }
    KRATOS_EXPECT_EQ(n_incised, 2);
    KRATOS_EXPECT_EQ(n_intersected, 2);

    // Check edge distances
    const auto &r_edge_dist_elem_2 = volume_part.GetElement(2).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_4 = volume_part.GetElement(4).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_elem_2 = {-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_3 = {-1.0,-1.0,0.5};
    const std::vector<double> expected_values_elem_4 = {0.5,0.5,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_2, expected_values_elem_2, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3, expected_values_elem_3, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4, expected_values_elem_4, 1.0e-6);

    //Check extra edge distances - elem_2 is not cut at all, elem_3 is incised, elem_4 is intersected
    const auto &r_edge_dist_elem_2_extra = volume_part.GetElement(2).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_4_extra = volume_part.GetElement(4).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_values_elem_2_extra = {-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_3_extra = {-1.0,0.5,-1.0};
    const std::vector<double> expected_values_elem_4_extra = {-1.0,-1.0,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_2_extra, expected_values_elem_2_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3_extra, expected_values_elem_3_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4_extra, expected_values_elem_4_extra, 1.0e-6);

    //Check elemental distances with extrapolated - values of elements that are not incised or intersected are characteristic length
    const auto &r_elem_dist_elem_2_extra = volume_part.GetElement(2).GetValue(ELEMENTAL_DISTANCES);
    const auto &r_elem_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    const auto &r_elem_dist_elem_4_extra = volume_part.GetElement(4).GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values_elem_2_extra_nodal = {1.41421,1.41421,1.41421};
    const std::vector<double> expected_values_elem_3_extra_nodal = {-0.25,0.25,0.25};
    const std::vector<double> expected_values_elem_4_extra_nodal = {-0.25,-0.25,0.25};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_2_extra, expected_values_elem_2_extra_nodal, 1.0e-5);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3_extra, expected_values_elem_3_extra_nodal, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_4_extra, expected_values_elem_4_extra_nodal, 1.0e-5);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessIncisedVsIntersected3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate tetrahedron elements
    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.CreateNewNode(1, -0.1, -0.1, -0.1);
    volume_part.CreateNewNode(2,  0.9, -0.1, -0.1);
    volume_part.CreateNewNode(3,  0.9,  0.9, -0.1);
    volume_part.CreateNewNode(4, -0.1,  0.9, -0.1);
    volume_part.CreateNewNode(5, -0.1, -0.1,  0.9);
    volume_part.CreateNewNode(6,  0.9, -0.1,  0.9);
    volume_part.CreateNewNode(7,  0.9,  0.9,  0.9);
    volume_part.CreateNewNode(8, -0.1,  0.9,  0.9);
    Properties::Pointer p_properties_0(new Properties(0));
    volume_part.CreateNewElement("Element3D4N", 1, {1, 6, 2, 4}, p_properties_0);
    volume_part.CreateNewElement("Element3D4N", 2, {1, 5, 6, 4}, p_properties_0);
    volume_part.CreateNewElement("Element3D4N", 3, {6, 5, 8, 4}, p_properties_0);
    volume_part.CreateNewElement("Element3D4N", 4, {2, 6, 3, 4}, p_properties_0);
    volume_part.CreateNewElement("Element3D4N", 5, {6, 7, 3, 4}, p_properties_0);
    volume_part.CreateNewElement("Element3D4N", 6, {6, 8, 7, 4}, p_properties_0);

    // Generate the skin line
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1,  0.3, -0.3, 1.0);
    skin_part.CreateNewNode(2,  0.3,  0.8, 1.0);
    skin_part.CreateNewNode(3, -0.3, -0.3, 0.0);
    skin_part.CreateNewNode(4, -0.3,  0.8, 0.0);
    Properties::Pointer p_properties_1(new Properties(1));
    skin_part.CreateNewElement("Element3D3N", 1, {{1,2,3}}, p_properties_1);
    skin_part.CreateNewElement("Element3D3N", 2, {{3,2,4}}, p_properties_1);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    parameters.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Count intersected and incised elements
    const uint8_t n_elements = 6;
    const uint8_t n_edges = 6;
    uint8_t n_cut_edges = 0;
    uint8_t n_incised = 0;
    uint8_t n_intersected = 0;
    for (uint8_t i_elem = 0; i_elem < n_elements; ++i_elem) {
        const auto &r_edge_dist = (volume_part.ElementsBegin() + i_elem)->GetValue(ELEMENTAL_EDGE_DISTANCES);
        n_cut_edges = 0;
        for (uint8_t j = 0; j < n_edges; ++j) {
            if (r_edge_dist[j] >= 0){
                n_cut_edges++;
            }
        }
        if (n_cut_edges > 0){
            if (n_cut_edges > 2){
                KRATOS_EXPECT_TRUE(i_elem + 1 == 2 || i_elem + 1 == 3);
                if (n_cut_edges == 3) {
                    bool is_incised = false;
                    for (std::size_t i = 0; i < n_edges; i++) {
                        double tolerance = std::numeric_limits<double>::epsilon();
                        const auto &r_edge_dist_extra = (volume_part.ElementsBegin() + i_elem)->GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
                        if ( std::abs(r_edge_dist_extra[i] - (-1.0)) > tolerance ) {
                            is_incised = true;
                        }
                    }
                    if (is_incised) {
                        n_incised++;
                    } else {
                        n_intersected++;
                    }
                } else {
                    n_intersected++;
                }
            } else {
                n_incised++;
                KRATOS_EXPECT_EQ(i_elem + 1, 6);
            }
        }
    }
    KRATOS_EXPECT_EQ(n_incised, 2);
    KRATOS_EXPECT_EQ(n_intersected, 1);

    // Check edge distances - elem_1,4,5 are not cut, elem_2 is intersected, elem_3 and elem_6 are (differently) incised
    const auto &r_edge_dist_elem_1 = volume_part.GetElement(1).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_2 = volume_part.GetElement(2).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_6 = volume_part.GetElement(6).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_elem_1 = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_2 = {0.433333,0.34,-1,-1,0.566667,-1};
    const std::vector<double> expected_values_elem_3 = {0.66,-1.0,0.34,-1.0,0.56666667,-1.0};
    const std::vector<double> expected_values_elem_6 = {0.66,-1.0,-1.0,-1.0,-1.0,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_1, expected_values_elem_1, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_2, expected_values_elem_2, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3, expected_values_elem_3, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_6, expected_values_elem_6, 1.0e-6);

    //Check extra edge distances - elem_2 and elem_3 both have three cut edges, but only elem_2 is completely intersected
    const auto &r_edge_dist_elem_1_extra = volume_part.GetElement(1).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_2_extra = volume_part.GetElement(2).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_6_extra = volume_part.GetElement(6).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_values_elem_1_extra = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_2_extra = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    const std::vector<double> expected_values_elem_3_extra = {-1,-1,-1,-1,-1,0.566667};
    const std::vector<double> expected_values_elem_6_extra = {-1,0.34,-1,-1,0.566667,-1};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_1_extra, expected_values_elem_1_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_2_extra, expected_values_elem_2_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3_extra, expected_values_elem_3_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_6_extra, expected_values_elem_6_extra, 1.0e-6);

    //Check elemental distances with extrapolated - values of elements that are not incised are characteristic length
    const auto &r_elem_dist_elem_2_extra = volume_part.GetElement(2).GetValue(ELEMENTAL_DISTANCES);
    const auto &r_elem_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    const auto &r_elem_dist_elem_6_extra = volume_part.GetElement(4).GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values_elem_2_extra_nodal = {-0.222948,0.291548,-0.565945,-0.222948};
    const std::vector<double> expected_values_elem_3_extra_nodal = {-0.565945,0.291548,0.291548,-0.222948};
    const std::vector<double> expected_values_elem_6_extra_nodal = {1.73205,1.73205,1.73205,1.73205};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_2_extra, expected_values_elem_2_extra_nodal, 1.0e-5);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3_extra, expected_values_elem_3_extra_nodal, 1.0e-5);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_6_extra, expected_values_elem_6_extra_nodal, 1.0e-5);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessEndAtEdge2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.25, 0.25, 0.0);
    skin_part.CreateNewNode(2,  0.5, 0.1, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    parameters.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Count intersected and incised elements
    const uint8_t n_elements = 8;
    const uint8_t n_edges = 3;
    uint8_t n_cut_edges = 0;
    uint8_t n_incised = 0;
    uint8_t n_intersected = 0;
    for (uint8_t i = 0; i < n_elements; ++i) {
        const auto &r_edge_dist = (volume_part.ElementsBegin() + i)->GetValue(ELEMENTAL_EDGE_DISTANCES);
        n_cut_edges = 0;
        for (uint8_t j = 0; j < n_edges; ++j) {
            if (r_edge_dist[j] >= 0){
                n_cut_edges++;
            }
        }
        if (n_cut_edges > 0){
            if (n_cut_edges > 1){
                n_intersected++;
            } else {
                n_incised++;
            }
        }
    }
    // Both Ends at Edges lead to an intersected element (elem_4 and elem_8), one leads to an incised element additionally (elem_3)
    KRATOS_EXPECT_EQ(n_incised, 1);
    KRATOS_EXPECT_EQ(n_intersected, 3);

    // Check edge distances - elements 1,2,5,6 are not cut at all
    const auto &r_edge_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_4 = volume_part.GetElement(4).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_7 = volume_part.GetElement(7).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_8 = volume_part.GetElement(8).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_elem_3 = {-1.0,-1.0,0.5};
    const std::vector<double> expected_values_elem_4 = {0.4,0.5,-1.0};
    const std::vector<double> expected_values_elem_7 = {-1.0,0.6,0.3333333};
    const std::vector<double> expected_values_elem_8 = {0.2,0.6666667,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3, expected_values_elem_3, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4, expected_values_elem_4, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_7, expected_values_elem_7, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_8, expected_values_elem_8, 1.0e-6);

    //Check extra edge distances - elem_3 is incised
    const auto &r_edge_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_values_elem_3_extra = {-1.0,0.4,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3_extra, expected_values_elem_3_extra, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessCutOnEdge2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    const double plane_height = 0.0;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.25, plane_height, 0.0);
    skin_part.CreateNewNode(2,  0.4, plane_height, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Count intersected and incised elements
    const uint8_t n_elements = 8;
    const uint8_t n_edges = 3;
    uint8_t n_cut_edges = 0;
    uint8_t n_incised = 0;
    uint8_t n_intersected = 0;
    for (uint8_t i = 0; i < n_elements; ++i) {
        const auto &r_edge_dist = (volume_part.ElementsBegin() + i)->GetValue(ELEMENTAL_EDGE_DISTANCES);
        n_cut_edges = 0;
        for (uint8_t j = 0; j < n_edges; ++j) {
            if (r_edge_dist[j] >= 0){
                n_cut_edges++;
            }
        }
        if (n_cut_edges > 0){
            if (n_cut_edges > 1){
                n_intersected++;
            } else {
                n_incised++;
            }
        }
    }

    KRATOS_EXPECT_EQ(n_incised, 2);
    KRATOS_EXPECT_EQ(n_intersected, 1);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessCutThroughNode2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.4,  0.2, 0.0);
    skin_part.CreateNewNode(2,  0.4, -0.2, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters;
    parameters.AddBool("calculate_elemental_edge_distances", true);
    parameters.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    parameters.AddBool("use_positive_epsilon_for_zero_values", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // Count intersected and incised elements
    const uint8_t n_elements = 8;
    const uint8_t n_edges = 3;
    uint8_t n_cut_edges = 0;
    uint8_t n_incised = 0;
    uint8_t n_intersected = 0;
    for (uint8_t i = 0; i < n_elements; ++i) {
        const auto &r_edge_dist = (volume_part.ElementsBegin() + i)->GetValue(ELEMENTAL_EDGE_DISTANCES);
        n_cut_edges = 0;
        for (uint8_t j = 0; j < n_edges; ++j) {
            if (r_edge_dist[j] >= 0){
                n_cut_edges++;
            }
        }
        if (n_cut_edges > 0){
            if (n_cut_edges > 1){
                n_intersected++;
            } else {
                n_incised++;
            }
        }
    }

    KRATOS_EXPECT_EQ(n_intersected, 4);
    KRATOS_EXPECT_EQ(n_incised, 2);

    // Check edge distances -> elem_4 and elem_5 are detected as only incised, elem_3 and elem_6 are incised
    const auto &r_edge_dist_elem_2 = volume_part.GetElement(2).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const auto &r_edge_dist_elem_4 = volume_part.GetElement(4).GetValue(ELEMENTAL_EDGE_DISTANCES);
    const std::vector<double> expected_values_elem_2 = {1.0,0.0,-1.0};
    const std::vector<double> expected_values_elem_3 = {-1.0,-1.0, 1.0 / 3.0};
    const std::vector<double> expected_values_elem_4 = {-1.0, 2.0 / 3.0 ,1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_2, expected_values_elem_2, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3, expected_values_elem_3, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4, expected_values_elem_4, 1.0e-6);

    //Check extra edge distances - elem_3 is incised, elem_4 is detected as intersected
    const auto &r_edge_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const auto &r_edge_dist_elem_4_extra = volume_part.GetElement(4).GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
    const std::vector<double> expected_values_elem_3_extra = {-1.0,0.5,-1.0};
    const std::vector<double> expected_values_elem_4_extra = {-1.0,-1.0,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3_extra, expected_values_elem_3_extra, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4_extra, expected_values_elem_4_extra, 1.0e-6);

    //Check elemental distances with extrapolated
    const auto &r_elem_dist_elem_3_extra = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    const auto &r_elem_dist_elem_4_extra = volume_part.GetElement(4).GetValue(ELEMENTAL_DISTANCES);
    const std::vector<double> expected_values_elem_3_extra_nodal = {-0.223607,0.447214,0.223607};
    const std::vector<double> expected_values_elem_4_extra_nodal = {-0.223607,0.0,0.447214};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3_extra, expected_values_elem_3_extra_nodal, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_4_extra, expected_values_elem_4_extra_nodal, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessFlags2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    const double plane_height = 0.25;
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.4, plane_height, 0.0);
    skin_part.CreateNewNode(2,  0.4, plane_height, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Both falgs are not given/ false
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc_0(volume_part, skin_part);
    disc_dist_proc_0.Execute();

    KRATOS_ERROR_IF(volume_part.ElementsBegin()->Has(ELEMENTAL_EDGE_DISTANCES));
    KRATOS_ERROR_IF(volume_part.ElementsBegin()->Has(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED));

    // Check elemental distances - elem_3 is incised
    auto& r_elem_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    std::vector<double> expected_values_elem_3 = {1.41421,1.41421,1.41421};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3, expected_values_elem_3, 1.0e-5);

    // Only CALCULATE_ELEMENTAL_EDGE_DISTANCES is given
    Parameters parameters_1;
    parameters_1.AddBool("calculate_elemental_edge_distances", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc_1(volume_part, skin_part, parameters_1);
    disc_dist_proc_1.Execute();

    KRATOS_ERROR_IF(volume_part.ElementsBegin()->Has(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED));

    // Check elemental distances - elem_3 is incised
    r_elem_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    expected_values_elem_3 = {1.41421,1.41421,1.41421};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3, expected_values_elem_3, 1.0e-5);

    // Check edge distances - elem_3 is incised
    auto &r_edge_dist = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    size_t n_cut_edges = 0;
    for (uint8_t j = 0; j < 3; ++j) {
        if (r_edge_dist[j] >= 0){
            n_cut_edges++;
        }
    }
    KRATOS_EXPECT_EQ(n_cut_edges, 1);

    //Both flags are given: CALCULATE_ELEMENTAL_EDGE_DISTANCES and CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED
    Parameters parameters_2;
    parameters_2.AddBool("calculate_elemental_edge_distances", true);
    parameters_2.AddBool("calculate_elemental_edge_distances_extrapolated", true);
    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc_2(volume_part, skin_part, parameters_2);
    disc_dist_proc_2.Execute();

    // Check elemental distances - elem_3 is incised
    r_elem_dist_elem_3 = volume_part.GetElement(3).GetValue(ELEMENTAL_DISTANCES);
    expected_values_elem_3 = {-0.25,0.25,0.25};
    KRATOS_EXPECT_VECTOR_NEAR(r_elem_dist_elem_3, expected_values_elem_3, 1.0e-5);

    // Check edge distances - elem_3 is incised
    r_edge_dist = volume_part.GetElement(3).GetValue(ELEMENTAL_EDGE_DISTANCES);
    n_cut_edges = 0;
    for (uint8_t j = 0; j < 3; ++j) {
        if (r_edge_dist[j] >= 0){
            n_cut_edges++;
        }
    }
    KRATOS_EXPECT_EQ(n_cut_edges, 1);
}

KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessCloseToVertexIntersection3D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the element
    ModelPart &fluid_part = current_model.CreateModelPart("Surface");
    fluid_part.CreateNewNode(1, 0.498262,    0.296646,    -0.0435666);
    fluid_part.CreateNewNode(2, 0.494408,    0.298003,    -0.0436762);
    fluid_part.CreateNewNode(3, 0.497984,    0.301717,    -0.046839);
    fluid_part.CreateNewNode(4, 0.496292,    0.295241,    -0.0478173);

    Properties::Pointer p_properties_0(new Properties(0));
    fluid_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(901, 0.2490485, -0.01, -0.02178895);
    skin_part.CreateNewNode(903, 0.498097, 2.0, -0.0435779);
    skin_part.CreateNewNode(904, 0.498097, -0.01, -0.0435779);
    skin_part.CreateNewNode(905, 2.0, 2.0, 0.0);

    Properties::Pointer p_properties = skin_part.CreateNewProperties(0);
    skin_part.CreateNewElement("Element3D3N", 901, { 901,904,903}, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 904,905,903}, p_properties);

    // Compute the discontinuous distance function
    CalculateDiscontinuousDistanceToSkinProcess<3> disc_dist_proc(fluid_part, skin_part);
    disc_dist_proc.Execute();
    auto p_elem = fluid_part.ElementsBegin();
    KRATOS_EXPECT_TRUE(p_elem->Is(TO_SPLIT));
}


KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceNewVariablese2D, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a fluid mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N"
    })");

    ModelPart& volume_part = current_model.CreateModelPart("Volume");
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    // Generate the skin line
    ModelPart& skin_part = current_model.CreateModelPart("Skin");
    skin_part.CreateNewNode(1, -0.25, 0.25, 0.0);
    skin_part.CreateNewNode(2,  0.5, 0.1, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties);

    // Compute the discontinuous distance function (including edge distances)
    Parameters parameters = Parameters(R"(
    {
        "elemental_distances_variable"                          : "EXTERNAL_FORCES_VECTOR",
        "elemental_edge_distances_variable"                     : "INTERNAL_FORCES_VECTOR",
        "elemental_edge_distances_extrapolated_variable"        : "CONTACT_FORCES_VECTOR",
        "calculate_elemental_edge_distances"                    : true,
        "calculate_elemental_edge_distances_extrapolated"       : true
    })" );

    CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(volume_part, skin_part, parameters);
    disc_dist_proc.Execute();

    // checking new ELEMENTAL_DISTANCES
    KRATOS_EXPECT_TRUE(volume_part.GetElement(3).Has(EXTERNAL_FORCES_VECTOR))
    // checking we dont have the default
    KRATOS_EXPECT_TRUE(!volume_part.GetElement(3).Has(ELEMENTAL_DISTANCES))
    // checking new ELEMENTAL_EDGE_DISTANCES
    KRATOS_EXPECT_TRUE(volume_part.GetElement(3).Has(INTERNAL_FORCES_VECTOR))
    // checking we dont have the default
    KRATOS_EXPECT_TRUE(!volume_part.GetElement(3).Has(ELEMENTAL_EDGE_DISTANCES))
    // checking new ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED
    KRATOS_EXPECT_TRUE(volume_part.GetElement(3).Has(INTERNAL_FORCES_VECTOR))
    // checking we dont have the default
    KRATOS_EXPECT_TRUE(!volume_part.GetElement(3).Has(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED))

    // Check edge distances - elements 1,2,5,6 are not cut at all
    const auto &r_edge_dist_elem_3 = volume_part.GetElement(3).GetValue(INTERNAL_FORCES_VECTOR);
    const auto &r_edge_dist_elem_4 = volume_part.GetElement(4).GetValue(INTERNAL_FORCES_VECTOR);
    const auto &r_edge_dist_elem_7 = volume_part.GetElement(7).GetValue(INTERNAL_FORCES_VECTOR);
    const auto &r_edge_dist_elem_8 = volume_part.GetElement(8).GetValue(INTERNAL_FORCES_VECTOR);
    const std::vector<double> expected_values_elem_3 = {-1.0,-1.0,0.5};
    const std::vector<double> expected_values_elem_4 = {0.4,0.5,-1.0};
    const std::vector<double> expected_values_elem_7 = {-1.0,0.6,0.3333333};
    const std::vector<double> expected_values_elem_8 = {0.2,0.6666667,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3, expected_values_elem_3, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_4, expected_values_elem_4, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_7, expected_values_elem_7, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_8, expected_values_elem_8, 1.0e-6);

    //Check extra edge distances - elem_3 is incised
    const auto &r_edge_dist_elem_3_extra = volume_part.GetElement(3).GetValue(CONTACT_FORCES_VECTOR);
    const std::vector<double> expected_values_elem_3_extra = {-1.0,0.4,-1.0};
    KRATOS_EXPECT_VECTOR_NEAR(r_edge_dist_elem_3_extra, expected_values_elem_3_extra, 1.0e-6);
}

}  // namespace Kratos::Testing.
