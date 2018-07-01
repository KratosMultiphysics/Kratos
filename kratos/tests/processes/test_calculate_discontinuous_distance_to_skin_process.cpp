//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessCubeInCube, KratosCoreFastSuite)
    {
        Model current_model;

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, -0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2,  0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3,  0.5,  0.5, -0.5);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, -0.5,  0.5, -0.5);
        Node<3>::Pointer p_point_5 = Kratos::make_shared<Node<3>>(5, -0.5, -0.5,  0.5);
        Node<3>::Pointer p_point_6 = Kratos::make_shared<Node<3>>(6,  0.5, -0.5,  0.5);
        Node<3>::Pointer p_point_7 = Kratos::make_shared<Node<3>>(7,  0.5,  0.5,  0.5);
        Node<3>::Pointer p_point_8 = Kratos::make_shared<Node<3>>(8, -0.5,  0.5,  0.5);

        Hexahedra3D8<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions":   5,
            "element_name":     "Element3D4N"
        })");

        ModelPart& volume_part = current_model.CreateModelPart("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_dist_elem_1 = (volume_part.ElementsBegin() + 7)->GetValue(ELEMENTAL_DISTANCES);
        const auto &r_dist_elem_2 = (volume_part.ElementsEnd() - 7)->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_dist_elem_1[0], -0.15, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_1[1], 0.05, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_1[2], 0.05, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_1[3], -0.15, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_2[0], -0.05, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_2[1], 0.15, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_2[2], 0.15, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_elem_2[3], -0.05, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessSharpCornerInCube, KratosCoreFastSuite)
    {
        Model current_model;

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.0, 0.0, 0.0);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 1.0, 0.0, 0.0);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 1.0, 1.0, 0.0);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 0.0, 1.0, 0.0);
        Node<3>::Pointer p_point_5 = Kratos::make_shared<Node<3>>(5, 0.0, 0.0, 1.0);
        Node<3>::Pointer p_point_6 = Kratos::make_shared<Node<3>>(6, 1.0, 0.0, 1.0);
        Node<3>::Pointer p_point_7 = Kratos::make_shared<Node<3>>(7, 1.0, 1.0, 1.0);
        Node<3>::Pointer p_point_8 = Kratos::make_shared<Node<3>>(8, 0.0, 1.0, 1.0);

        Hexahedra3D8<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_dist_begin = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        const auto &r_dist_end = (volume_part.ElementsEnd() - 1)->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_dist_begin[0], 1.0, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_begin[1], 1.0, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_begin[2], 1.0, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_begin[3], 1.0, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_end[0], -0.406059, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_end[1], -0.489839, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_end[2], 0.388306, 1e-6);
        KRATOS_CHECK_NEAR(r_dist_end[3], 0.0649427, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessHorizontalPlane, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], 0.0714286, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[1], 0.0714286, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[2], 0.0714286, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[3], -0.0714286, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationSkewed, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], -0.569400, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[1], 0.1044875, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[2], -0.266495, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[3], 0.104487, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationVertical, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const Vector elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(elem_dist[0], -0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[1],  0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[2], -0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[3], -0.5, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessOneEdgeIntersection, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check elemental distance values
        const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], 1.0, 1e-10);
        KRATOS_CHECK_NEAR(r_elem_dist[1], 1.0, 1e-10);
        KRATOS_CHECK_NEAR(r_elem_dist[2], 1.0, 1e-10);
        KRATOS_CHECK_NEAR(r_elem_dist[3], 1.0, 1e-10);


    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessMultipleIntersections, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check elemental distance values
        const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], -0.0636738, 1e-7);
        KRATOS_CHECK_NEAR(r_elem_dist[1],  0.0342287, 1e-7);
        KRATOS_CHECK_NEAR(r_elem_dist[2], -0.0709816, 1e-7);
        KRATOS_CHECK_NEAR(r_elem_dist[3], -0.0159295, 1e-7);

    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessStandard, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], -0.108523, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[1], 0.0485713, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[2], 0.0764035, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[3], 0.0222594, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessBoundaryIntersection, KratosCoreFastSuite)
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

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const auto &r_elem_dist = volume_part.ElementsBegin()->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], -0.0984855, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[1], -0.00883326, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[2], 0.0352186, 1e-6);
        KRATOS_CHECK_NEAR(r_elem_dist[3], -0.103167, 1e-6);
    }

}  // namespace Testing.
}  // namespace Kratos.
