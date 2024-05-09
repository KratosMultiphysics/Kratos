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
#include "containers/model.h"
#include "includes/expect.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/apply_ray_casting_process.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(RayCastingProcessQuadrilateral2D, KratosCoreFastSuite)
{
    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.00, 0.00, 0.00);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.00, 10.00, 0.00);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 10.00, 10.00, 0.00);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 10.00, 0.00, 0.00);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   7,
        "element_name":     "Element2D3N"
    })");

    Model current_model;
    ModelPart &surface_part = current_model.CreateModelPart("Volume");
    surface_part.AddNodalSolutionStepVariable(VELOCITY);
    surface_part.AddNodalSolutionStepVariable(DISTANCE);
    surface_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
    StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

    for(auto& node : surface_part.Nodes()){
        node.GetSolutionStepValue(DISTANCE) = 1.00;
    }

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(901, 2.4, 3.4, 0.0);
    skin_part.CreateNewNode(902, 7.6, 3.4, 0.0);
    skin_part.CreateNewNode(903, 7.6, 6.6, 0.0);
    skin_part.CreateNewNode(904, 2.4, 6.6, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 901, {{901,902}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 902, {{902,903}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 903, {{903,904}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 904, {{904,901}}, p_properties);

    // Compute distance
    ApplyRayCastingProcess<2>(surface_part, skin_part).Execute();

    // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/surface_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_fluid.InitializeMesh(0.00);
    // gid_io_fluid.WriteMesh(surface_part.GetMesh());
    // gid_io_fluid.FinalizeMesh();
    // gid_io_fluid.InitializeResults(0, surface_part.GetMesh());
    // gid_io_fluid.WriteNodalResults(DISTANCE, surface_part.Nodes(), 0, 0);
    // gid_io_fluid.FinalizeResults();

    // GidIO<> gid_io_skin("/home/rzorrilla/Desktop/skin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_skin.InitializeMesh(0.00);
    // gid_io_skin.WriteMesh(skin_part.GetMesh());
    // gid_io_skin.FinalizeMesh();
    // gid_io_skin.InitializeResults(0, skin_part.GetMesh());
    // gid_io_skin.FinalizeResults();

    KRATOS_EXPECT_NEAR((surface_part.pGetNode(21))->FastGetSolutionStepValue(DISTANCE), -1.00, 1e-6);
    KRATOS_EXPECT_NEAR((surface_part.pGetNode(22))->FastGetSolutionStepValue(DISTANCE),  1.00, 1e-6);
    KRATOS_EXPECT_NEAR((surface_part.pGetNode(30))->FastGetSolutionStepValue(DISTANCE),  1.00, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RayCastingProcessSquareRing2D, KratosCoreFastSuite)
{
    // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.00, 0.00, 0.00);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.00, 10.00, 0.00);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 10.00, 10.00, 0.00);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 10.00, 0.00, 0.00);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   12,
        "element_name":     "Element2D3N"
    })");

    Model current_model;
    ModelPart &surface_part = current_model.CreateModelPart("Volume");
    surface_part.AddNodalSolutionStepVariable(VELOCITY);
    surface_part.AddNodalSolutionStepVariable(DISTANCE);
    surface_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
    StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

    for(auto& node : surface_part.Nodes()){
        node.GetSolutionStepValue(DISTANCE) = 1.00;
    }

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(901, 2.4, 2.4, 0.0);
    skin_part.CreateNewNode(902, 7.6, 2.4, 0.0);
    skin_part.CreateNewNode(903, 7.6, 7.6, 0.0);
    skin_part.CreateNewNode(904, 2.4, 7.6, 0.0);
    skin_part.CreateNewNode(905, 3.9, 3.9, 0.0);
    skin_part.CreateNewNode(906, 6.1, 3.9, 0.0);
    skin_part.CreateNewNode(907, 6.1, 6.1, 0.0);
    skin_part.CreateNewNode(908, 3.9, 6.1, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element2D2N", 901, {{901,902}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 902, {{902,903}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 903, {{903,904}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 904, {{904,901}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 905, {{905,906}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 906, {{906,907}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 907, {{907,908}}, p_properties);
    skin_part.CreateNewElement("Element2D2N", 908, {{908,905}}, p_properties);

    // Compute distance
    ApplyRayCastingProcess<2>(surface_part, skin_part).Execute();

    KRATOS_EXPECT_NEAR((surface_part.pGetNode(86))->FastGetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR((surface_part.pGetNode(88))->FastGetSolutionStepValue(DISTANCE), -1.00, 1e-6);
    KRATOS_EXPECT_NEAR((surface_part.pGetNode(112))->FastGetSolutionStepValue(DISTANCE), -1.00, 1e-6);
    KRATOS_EXPECT_NEAR((surface_part.pGetNode(138))->FastGetSolutionStepValue(DISTANCE),  1.00, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(TetrahedraInCubeRayCastingProcess, KratosCoreFastSuite)
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

    Hexahedra3D8<Node> geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":   10,
        "element_name":     "Element3D4N"
    })");

    Model current_model;
    ModelPart &volume_part = current_model.CreateModelPart("Volume");
    volume_part.AddNodalSolutionStepVariable(VELOCITY);
    volume_part.AddNodalSolutionStepVariable(DISTANCE);
    volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
    StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

    for(auto& node : volume_part.Nodes()){
        node.GetSolutionStepValue(DISTANCE) = -1.00;
    }

    // Generate the skin
    ModelPart &skin_part = current_model.CreateModelPart("Skin");
    skin_part.AddNodalSolutionStepVariable(VELOCITY);
    skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
    skin_part.CreateNewNode(902, 6.0, 2.0, 2.0);
    skin_part.CreateNewNode(903, 4.0, 6.0, 2.0);
    skin_part.CreateNewNode(904, 4.0, 4.0, 7.0);
    Properties::Pointer p_properties(new Properties(0));
    skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
    skin_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);

    // Compute distance
    ApplyRayCastingProcess<3>(volume_part, skin_part).Execute();


    Tetrahedra3D4<Node> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

    Point dummy(0.0,0.0,0.0);
    for(auto& node : volume_part.Nodes()){
        if(tetrahedra.IsInside(node.Coordinates(),dummy)){
            KRATOS_EXPECT_NEAR(node.GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
        }
    }

    // Note that we cannot check the outside because on the interface is not well defined
    KRATOS_EXPECT_NEAR(volume_part.GetNode(135).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetNode(136).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetNode(137).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetNode(256).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetNode(257).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
    KRATOS_EXPECT_NEAR(volume_part.GetNode(258).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);


    //GidIO<> gid_io_fluid("C:/Temp/Tests/distance_test_fluid", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    //gid_io_fluid.InitializeMesh(0.00);
    //gid_io_fluid.WriteMesh(volume_part.GetMesh());
    //gid_io_fluid.FinalizeMesh();
    //gid_io_fluid.InitializeResults(0, volume_part.GetMesh());
    //gid_io_fluid.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
    //gid_io_fluid.FinalizeResults();

}

}  // namespace Kratos::Python.
