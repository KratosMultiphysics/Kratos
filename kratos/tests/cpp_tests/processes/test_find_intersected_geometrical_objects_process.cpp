//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
// #include "includes/gid_io.h"
// #include "input_output/vtk_output.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/find_intersected_geometrical_objects_with_obb_process.h"
#include "processes/skin_detection_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcess2D, KratosCoreFastSuite)
        {
            NodeType::Pointer p_point1(new NodeType(1, 0.0, 0.0, 0.0));
            NodeType::Pointer p_point2(new NodeType(2, 0.0, 1.0, 0.0));
            NodeType::Pointer p_point3(new NodeType(3, 1.0, 1.0, 0.0));
            NodeType::Pointer p_point4(new NodeType(4, 1.0, 0.0, 0.0));

            Quadrilateral2D4<NodeType > geometry(p_point1, p_point2, p_point3, p_point4);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions" : 3,
                "element_name" : "Element2D3N"
            }  )");

            Model current_model;
            ModelPart &surface_part = current_model.CreateModelPart("Surface");
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            skin_part.CreateNewNode(100, -0.3, 0.5, 0.0);
            skin_part.CreateNewNode(200, 0.6, 0.5, 0.0);
            Properties::Pointer p_properties(new Properties(0));
            skin_part.CreateNewElement("Element2D2N", 1, {{ 100,200 }}, p_properties);
            StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();
            FindIntersectedGeometricalObjectsProcess find_intersections(surface_part, skin_part);
            find_intersections.Execute();

            // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/surface_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io_fluid.InitializeMesh(0.00);
            // gid_io_fluid.WriteMesh(surface_part.GetMesh());
            // gid_io_fluid.FinalizeMesh();
            // gid_io_fluid.InitializeResults(0, surface_part.GetMesh());
            // gid_io_fluid.FinalizeResults();

            // GidIO<> gid_io_skin("/home/rzorrilla/Desktop/skin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io_skin.InitializeMesh(0.00);
            // gid_io_skin.WriteMesh(skin_part.GetMesh());
            // gid_io_skin.FinalizeMesh();
            // gid_io_skin.InitializeResults(0, skin_part.GetMesh());
            // gid_io_skin.FinalizeResults();

            KRATOS_CHECK((surface_part.Elements()[3]).Is(SELECTED));
            KRATOS_CHECK((surface_part.Elements()[4]).Is(SELECTED));
            KRATOS_CHECK((surface_part.Elements()[9]).Is(SELECTED));
            KRATOS_CHECK((surface_part.Elements()[10]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[0]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[1]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[2]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[5]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[6]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[7]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[8]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[11]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[12]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[13]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[14]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[15]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[16]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[17]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcess2D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");
            r_surface_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_surface_part.CreateNewNode(2, 0.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(4, 1.0, 0.0, 0.0);
            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");
            r_skin_part.CreateNewNode(5, -0.3, 0.5, 0.0);
            r_skin_part.CreateNewNode(6, 0.6, 0.5, 0.0);
            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_properties_0);
            r_skin_part.CreateNewCondition("LineCondition2D2N", 5, {{ 5,6 }}, p_properties_1);
            FindIntersectedGeometricalObjectsProcess find_intersections(r_surface_part, r_skin_part);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcessNoIntersection2D, KratosCoreFastSuite)
        {
            NodeType::Pointer p_point1(new NodeType(1, 0.0, 0.0, 0.0));
            NodeType::Pointer p_point2(new NodeType(2, 0.0, 1.0, 0.0));
            NodeType::Pointer p_point3(new NodeType(3, 1.0, 1.0, 0.0));
            NodeType::Pointer p_point4(new NodeType(4, 1.0, 0.0, 0.0));

            Quadrilateral2D4<NodeType > geometry(p_point1, p_point2, p_point3, p_point4);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions" : 3,
                "element_name" : "Element2D3N"
            }  )");

              Model current_model;
            ModelPart &surface_part = current_model.CreateModelPart("Surface");
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            skin_part.CreateNewNode(100, 0.3, -0.5, 0.0);
            skin_part.CreateNewNode(200, 0.6, -0.5, 0.0);
            Properties::Pointer p_properties(new Properties(0));
            skin_part.CreateNewElement("Element2D2N", 1, {{ 100,200 }}, p_properties);
            StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();
            FindIntersectedGeometricalObjectsProcess find_intersections(surface_part, skin_part);
            find_intersections.Execute();

            for (auto it_elem = surface_part.ElementsBegin(); it_elem != surface_part.ElementsEnd(); ++it_elem){
                KRATOS_CHECK_IS_FALSE(it_elem->Is(SELECTED));
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessNoIntersection2D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");
            r_surface_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_surface_part.CreateNewNode(2, 0.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(4, 1.0, 0.0, 0.0);
            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");
            r_skin_part.CreateNewNode(5, -0.1, 0.0, 0.0);
            r_skin_part.CreateNewNode(6, -0.1, 1.0, 0.0);
            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_properties_0);
            r_skin_part.CreateNewCondition("LineCondition2D2N", 5, {{ 5,6 }}, p_properties_1);
            FindIntersectedGeometricalObjectsProcess find_intersections(r_surface_part, r_skin_part);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[3]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[4]).IsNot(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessIntersectionExtendedOBB2D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");
            r_surface_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_surface_part.CreateNewNode(2, 0.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_surface_part.CreateNewNode(4, 1.0, 0.0, 0.0);
            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");
            r_skin_part.CreateNewNode(5, -0.1, 0.0, 0.0);
            r_skin_part.CreateNewNode(6, -0.1, 1.0, 0.0);
            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_properties_0);
            r_surface_part.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_properties_0);
            r_skin_part.CreateNewCondition("LineCondition2D2N", 5, {{ 5,6 }}, p_properties_1);

            Parameters parameters = Parameters(R"(
            {
                "intersected_model_part_name"  : "Main.Surface",
                "intersecting_model_part_name" : "Main.Boundaries",
                "bounding_box_factor"          : 0.2,
                "debug_obb"                    : false
            })" );

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(current_model, parameters);
//             FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.2);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[3]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[4]).IsNot(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessIntersectionExtendedOBB3D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            r_surface_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            r_surface_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
            r_surface_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);

            // Now we create the "elements"

            r_surface_part.CreateNewElement("Element3D4N", 1, {{5, 3, 8, 6}}, p_properties_0);
            r_surface_part.CreateNewElement("Element3D4N", 2, {{4, 6, 7, 3}}, p_properties_0);
            r_surface_part.CreateNewElement("Element3D4N", 3, {{2, 3, 5, 6}}, p_properties_0);
            r_surface_part.CreateNewElement("Element3D4N", 4, {{7, 8, 3, 6}}, p_properties_0);
            r_surface_part.CreateNewElement("Element3D4N", 5, {{4, 1, 6, 3}}, p_properties_0);
            r_surface_part.CreateNewElement("Element3D4N", 6, {{3, 2, 1, 6}}, p_properties_0);

            // Generate skin
            Parameters surface_parameters = Parameters(R"(
            {
                "name_auxiliar_model_part"              : "Main.Surface",
                "name_auxiliar_condition"               : "Condition",
                "list_model_parts_to_assign_conditions" : [],
                "echo_level"                            : 0
            })");

//             SkinDetectionProcess<3> surface_skin_process(r_surface_part, surface_parameters);
//             surface_skin_process.Execute();

            // Deterministic creation
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",1, {{ 1, 6, 2 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",2, {{ 3, 1, 2 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",3, {{ 5, 6, 8 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",4, {{ 4, 7, 6 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",5, {{ 5, 8, 3 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",6, {{ 2, 5, 3 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",7, {{ 7, 3, 8 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",8, {{ 2, 6, 5 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",9, {{ 7, 8, 6 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",10, {{ 4, 6, 1 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",11, {{ 4, 3, 7 }}, p_properties_1);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N",12, {{ 4, 1, 3 }}, p_properties_1);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");
            r_skin_part.CreateNewNode(9 , -0.1 , 1.0 , 1.0);
            r_skin_part.CreateNewNode(10 , -0.1 , 1.0 , 0.0);
            r_skin_part.CreateNewNode(11 , -0.1 , 0.0 , 1.0);
            r_skin_part.CreateNewNode(12 , -0.1 , 0.0 , 0.0);

            const std::size_t number_of_conditions = r_main_model_part.NumberOfConditions();
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", number_of_conditions + 1, {{ 9, 10, 11 }}, p_properties_1);
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", number_of_conditions + 2, {{ 10, 11, 12 }}, p_properties_1);

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.1);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[3]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[4]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[5]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[6]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[7]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[8]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[9]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[10]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[11]).IsNot(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[12]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessTriangleTriangleOBB3D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{3, 4, 1}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 3, {{ 5, 6, 7 }}, p_properties_1);
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 4, {{ 7, 8, 5 }}, p_properties_1);

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.1);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessTriangleTriangleOBB3DOrthogonalBase, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{3, 4, 1}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 3, {{ 5, 6, 7 }}, p_properties_1);
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 4, {{ 7, 8, 5 }}, p_properties_1);

            Parameters intersect_parameters = Parameters(R"(
            {
                "intersected_model_part_name"  : "Main.Surface",
                "intersecting_model_part_name" : "Main.Boundaries",
                "bounding_box_factor"          : 0.1,
                "debug_obb"                    : false,
                "OBB_intersection_type"        : "SeparatingAxisTheorem",
                "build_from_bounding_box"      : false,
                "intersecting_conditions"      : true,
                "intersecting_elements"        : true,
                "intersected_conditions"       : true,
                "intersected_elements"         : true
            })" );

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(current_model, intersect_parameters);
            find_intersections.Execute();

//             Parameters vtk_parameters = Parameters(R"(
//             {
//                 "model_part_name"                    : "Main",
//                 "file_format"                        : "ascii",
//                 "output_precision"                   : 7,
//                 "output_control_type"                : "step",
//                 "output_frequency"                   : 1.0,
//                 "output_sub_model_parts"             : true,
//                 "folder_name"                        : "VTK_Output",
//                 "custom_name_prefix"                 : "",
//                 "save_output_files_in_folder"        : false,
//                 "write_deformed_configuration"       : false,
//                 "write_properties_id"                : false,
//                 "nodal_solution_step_data_variables" : [],
//                 "nodal_data_value_variables"         : [],
//                 "nodal_flags"                        : [],
//                 "element_data_value_variables"       : [],
//                 "element_flags"                      : [],
//                 "condition_data_value_variables"     : [],
//                 "condition_flags"                    : ["SELECTED"],
//                 "gauss_point_variables"              : []
//             })" );
//
//             VtkOutput debug(r_main_model_part, vtk_parameters);
//             debug.PrintOutput();
//
//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessQuadrilateralQuadrilateralOBB3D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D4N", 1, {{1, 2, 3, 4}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D4N", 2, {{ 5, 6, 7, 8 }}, p_properties_1);

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.1);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessQuadrilateralQuadrilateralOBB3DOrthogonalBase, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D4N", 1, {{1, 2, 3, 4}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D4N", 2, {{ 5, 6, 7, 8 }}, p_properties_1);

            Parameters intersect_parameters = Parameters(R"(
            {
                "intersected_model_part_name"  : "Main.Surface",
                "intersecting_model_part_name" : "Main.Boundaries",
                "bounding_box_factor"          : 0.1,
                "debug_obb"                    : false,
                "OBB_intersection_type"        : "SeparatingAxisTheorem",
                "build_from_bounding_box"      : false,
                "intersecting_conditions"      : true,
                "intersecting_elements"        : true,
                "intersected_conditions"       : true,
                "intersected_elements"         : true
            })" );

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(current_model, intersect_parameters);
            find_intersections.Execute();

//             Parameters vtk_parameters = Parameters(R"(
//             {
//                 "model_part_name"                    : "Main",
//                 "file_format"                        : "ascii",
//                 "output_precision"                   : 7,
//                 "output_control_type"                : "step",
//                 "output_frequency"                   : 1.0,
//                 "output_sub_model_parts"             : true,
//                 "folder_name"                        : "VTK_Output",
//                 "custom_name_prefix"                 : "",
//                 "save_output_files_in_folder"        : false,
//                 "write_deformed_configuration"       : false,
//                 "write_properties_id"                : false,
//                 "nodal_solution_step_data_variables" : [],
//                 "nodal_data_value_variables"         : [],
//                 "nodal_flags"                        : [],
//                 "element_data_value_variables"       : [],
//                 "element_flags"                      : [],
//                 "condition_data_value_variables"     : [],
//                 "condition_flags"                    : ["SELECTED"],
//                 "gauss_point_variables"              : []
//             })" );
//
//             VtkOutput debug(r_main_model_part, vtk_parameters);
//             debug.PrintOutput();
//
//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessTriangleQuadrilateralOBB3D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D4N", 1, {{1, 2, 3, 4}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{ 5, 6, 7 }}, p_properties_1);
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 3, {{ 7, 8, 5 }}, p_properties_1);

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.1);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
        }


        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessTriangleQuadrilateralOBB3DOrthogonalBase, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D4N", 1, {{1, 2, 3, 4}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{ 5, 6, 7 }}, p_properties_1);
            r_skin_part.CreateNewCondition("SurfaceCondition3D3N", 3, {{ 7, 8, 5 }}, p_properties_1);

            Parameters intersect_parameters = Parameters(R"(
            {
                "intersected_model_part_name"  : "Main.Surface",
                "intersecting_model_part_name" : "Main.Boundaries",
                "bounding_box_factor"          : 0.1,
                "debug_obb"                    : false,
                "OBB_intersection_type"        : "SeparatingAxisTheorem",
                "build_from_bounding_box"      : false,
                "intersecting_conditions"      : true,
                "intersecting_elements"        : true,
                "intersected_conditions"       : true,
                "intersected_elements"         : true
            })" );

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(current_model, intersect_parameters);
            find_intersections.Execute();

//             Parameters vtk_parameters = Parameters(R"(
//             {
//                 "model_part_name"                    : "Main",
//                 "file_format"                        : "ascii",
//                 "output_precision"                   : 7,
//                 "output_control_type"                : "step",
//                 "output_frequency"                   : 1.0,
//                 "output_sub_model_parts"             : true,
//                 "folder_name"                        : "VTK_Output",
//                 "custom_name_prefix"                 : "",
//                 "save_output_files_in_folder"        : false,
//                 "write_deformed_configuration"       : false,
//                 "write_properties_id"                : false,
//                 "nodal_solution_step_data_variables" : [],
//                 "nodal_data_value_variables"         : [],
//                 "nodal_flags"                        : [],
//                 "element_data_value_variables"       : [],
//                 "element_flags"                      : [],
//                 "condition_data_value_variables"     : [],
//                 "condition_flags"                    : ["SELECTED"],
//                 "gauss_point_variables"              : []
//             })" );
//
//             VtkOutput debug(r_main_model_part, vtk_parameters);
//             debug.PrintOutput();
//
//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessQuadrilateralTriangleOBB3D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{3, 4, 1}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D4N", 3, {{ 5, 6, 7, 8 }}, p_properties_1);

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(r_surface_part, r_skin_part, 0.1);
            find_intersections.Execute();

//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedConditionsProcessQuadrilateralTriangleOBB3DOrthogonalBase, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_main_model_part = current_model.CreateModelPart("Main");
            r_main_model_part.GetProcessInfo().SetValue(STEP, 1);
            r_main_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
            ModelPart& r_surface_part = r_main_model_part.CreateSubModelPart("Surface");

            Properties::Pointer p_properties_0 = Kratos::make_shared<Properties>(0);
            Properties::Pointer p_properties_1 = Kratos::make_shared<Properties>(1);

            // Generate the nodes of the surface
            r_surface_part.CreateNewNode(1 , 0.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(2 , 1.0 , 0.0 , 1.0);
            r_surface_part.CreateNewNode(3 , 1.0 , 1.0 , 1.0);
            r_surface_part.CreateNewNode(4 , 0.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_properties_0);
            r_surface_part.CreateNewCondition("SurfaceCondition3D3N", 2, {{3, 4, 1}}, p_properties_0);

            ModelPart& r_skin_part = r_main_model_part.CreateSubModelPart("Boundaries");

            // Generate the nodes of the skin
            r_skin_part.CreateNewNode(5 , 0.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(6 , 1.0 , 0.0 , 1.01);
            r_skin_part.CreateNewNode(7 , 1.0 , 1.0 , 1.01);
            r_skin_part.CreateNewNode(8 , 0.0 , 1.0 , 1.01);

            // Now we create the "conditions"
            r_skin_part.CreateNewCondition("SurfaceCondition3D4N", 3, {{ 5, 6, 7, 8 }}, p_properties_1);

            Parameters intersect_parameters = Parameters(R"(
            {
                "intersected_model_part_name"  : "Main.Surface",
                "intersecting_model_part_name" : "Main.Boundaries",
                "bounding_box_factor"          : 0.1,
                "debug_obb"                    : false,
                "OBB_intersection_type"        : "SeparatingAxisTheorem",
                "build_from_bounding_box"      : false,
                "intersecting_conditions"      : true,
                "intersecting_elements"        : true,
                "intersected_conditions"       : true,
                "intersected_elements"         : true
            })" );

            FindIntersectedGeometricalObjectsWithOBBProcess find_intersections(current_model, intersect_parameters);
            find_intersections.Execute();

//             Parameters vtk_parameters = Parameters(R"(
//             {
//                 "model_part_name"                    : "Main",
//                 "file_format"                        : "ascii",
//                 "output_precision"                   : 7,
//                 "output_control_type"                : "step",
//                 "output_frequency"                   : 1.0,
//                 "output_sub_model_parts"             : true,
//                 "folder_name"                        : "VTK_Output",
//                 "custom_name_prefix"                 : "",
//                 "save_output_files_in_folder"        : false,
//                 "write_deformed_configuration"       : false,
//                 "write_properties_id"                : false,
//                 "nodal_solution_step_data_variables" : [],
//                 "nodal_data_value_variables"         : [],
//                 "nodal_flags"                        : [],
//                 "element_data_value_variables"       : [],
//                 "element_flags"                      : [],
//                 "condition_data_value_variables"     : [],
//                 "condition_flags"                    : ["SELECTED"],
//                 "gauss_point_variables"              : []
//             })" );
//
//             VtkOutput debug(r_main_model_part, vtk_parameters);
//             debug.PrintOutput();
//
//             GidIO<> gid_io("test", GiD_PostBinary, SingleFile, WriteDeformed, WriteConditions);
//             gid_io.InitializeMesh(0.0);
//             gid_io.WriteMesh(r_main_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(0, r_main_model_part.GetMesh());
//             gid_io.FinalizeResults();

            KRATOS_CHECK((r_surface_part.Conditions()[1]).Is(SELECTED));
            KRATOS_CHECK((r_surface_part.Conditions()[2]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcessBoundingBoxIntersection2D, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            ModelPart &surface_part = current_model.CreateModelPart("Surface");

            // Create surface part
            NodeType::Pointer p_point1(new NodeType(1, 0.0, 0.0, 0.0));
            NodeType::Pointer p_point2(new NodeType(2, 0.0, 1.0, 0.0));
            NodeType::Pointer p_point3(new NodeType(3, 1.0, 1.0, 0.0));
            NodeType::Pointer p_point4(new NodeType(4, 1.0, 0.0, 0.0));
            Quadrilateral2D4<NodeType > geometry(p_point1, p_point2, p_point3, p_point4);
            Parameters mesher_parameters(R"({
                "number_of_divisions" : 3,
                "element_name" : "Element2D3N"
            })");
            StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

            // Create skin part
            skin_part.CreateNewNode(100, 0.2, 0.0, 0.0);
            skin_part.CreateNewNode(200, 0.5, 0.0, 0.0);
            Properties::Pointer p_properties(new Properties(0));
            skin_part.CreateNewElement("Element2D2N", 1, {{ 100,200 }}, p_properties);

            // Create and call the FindIntersectedGeometricalObjectsProcess
            FindIntersectedGeometricalObjectsProcess find_intersections(surface_part, skin_part);
            find_intersections.Execute();

            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[0]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[5]).Is(SELECTED));
            KRATOS_CHECK_IS_FALSE((surface_part.Elements()[6]).Is(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcess3D, KratosCoreFastSuite)
        {
            NodeType::Pointer p_point1(new NodeType(1, 0.00, 0.00, 0.00));
            NodeType::Pointer p_point2(new NodeType(2, 10.00, 0.00, 0.00));
            NodeType::Pointer p_point3(new NodeType(3, 10.00, 10.00, 0.00));
            NodeType::Pointer p_point4(new NodeType(4, 0.00, 10.00, 0.00));
            NodeType::Pointer p_point5(new NodeType(5, 0.00, 0.00, 10.00));
            NodeType::Pointer p_point6(new NodeType(6, 10.00, 0.00, 10.00));
            NodeType::Pointer p_point7(new NodeType(7, 10.00, 10.00, 10.00));
            NodeType::Pointer p_point8(new NodeType(8, 0.00, 10.00, 10.00));

            Hexahedra3D8<NodeType > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":2,
                "element_name": "Element3D4N"
            }  )");

              Model current_model;
            ModelPart &volume_part = current_model.CreateModelPart("Volume");
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            skin_part.CreateNewNode(1, 1., .2, 0.);
            skin_part.CreateNewNode(2, 1., .1, .5);
            skin_part.CreateNewNode(3, 1., .1, 0.);
            Properties::Pointer p_properties(new Properties(0));
            skin_part.CreateNewElement("Element3D3N", 1, { 1,2,3 }, p_properties);
            StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();
            FindIntersectedGeometricalObjectsProcess(volume_part, skin_part).Execute();
            KRATOS_CHECK(volume_part.GetElement(3).IsNot(SELECTED));
            KRATOS_CHECK(volume_part.GetElement(4).IsNot(SELECTED));
            KRATOS_CHECK(volume_part.GetElement(5).Is(SELECTED));
            KRATOS_CHECK(volume_part.GetElement(6).Is(SELECTED));
            for (std::size_t i = 7; i < volume_part.NumberOfElements(); i++)
                KRATOS_CHECK(volume_part.GetElement(i).IsNot(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcessNoIntersection3D, KratosCoreFastSuite)
        {
            // Generate the tetrahedron element
            Model current_model;
            ModelPart &volume_part = current_model.CreateModelPart("Volume");
            volume_part.CreateNewNode(34, 0.865646, 0.657938, 0.222985);
            volume_part.CreateNewNode(58, 0.770744, 0.570027, 0.204129);
            volume_part.CreateNewNode(73, 0.860052, 0.477371, 0.22713);
            volume_part.CreateNewNode(96, 0.803174, 0.485159, 0.326767);
            Properties::Pointer p_properties_0(new Properties(0));
            volume_part.CreateNewElement("Element3D4N", 139, {34, 58, 73, 96}, p_properties_0);

            // Generate the skin model part
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            skin_part.CreateNewNode(662, 0.766593, 0.532174, 0.275516);
            skin_part.CreateNewNode(723, 0.793214, 0.506089, 0.308981);
            skin_part.CreateNewNode(737, 0.794158, 0.544627, 0.315665);
            skin_part.CreateNewNode(801, 0.81563, 0.518347, 0.349863);
            skin_part.CreateNewNode(814, 0.811818, 0.567485, 0.356072);
            skin_part.CreateNewNode(777, 0.809491, 0.469669, 0.339392);
            skin_part.CreateNewNode(710, 0.7901, 0.455512, 0.309309);
            skin_part.CreateNewNode(682, 0.768283, 0.578834, 0.289503);
            skin_part.CreateNewNode(741, 0.786372, 0.593624, 0.321883);
            skin_part.CreateNewNode(652, 0.766584, 0.482207, 0.273911);
            Properties::Pointer p_properties_1(new Properties(1));
            skin_part.CreateNewElement("Element3D3N", 477, {662,723,737}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 478, {737,723,801}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 479, {737,801,814}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 480, {801,723,777}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 510, {777,723,710}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 467, {682,662,737}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 484, {737,741,682}, p_properties_1);
            skin_part.CreateNewElement("Element3D3N", 496, {723,652,710}, p_properties_1);

            // Call the intersections process
            FindIntersectedGeometricalObjectsProcess(volume_part, skin_part).Execute();

            // Check that there is no intersection
            KRATOS_CHECK(volume_part.GetElement(139).IsNot(SELECTED));
        }

        KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsConditionsProcessNoIntersection3D, KratosCoreFastSuite)
        {
            // Generate the tetrahedron element
            Model current_model;
            ModelPart &volume_part = current_model.CreateModelPart("Volume");
            volume_part.CreateNewNode(34, 0.865646, 0.657938, 0.222985);
            volume_part.CreateNewNode(58, 0.770744, 0.570027, 0.204129);
            volume_part.CreateNewNode(73, 0.860052, 0.477371, 0.22713);
            volume_part.CreateNewNode(96, 0.803174, 0.485159, 0.326767);
            Properties::Pointer p_properties_0(new Properties(0));
            volume_part.CreateNewElement("Element3D4N", 139, {34, 58, 73, 96}, p_properties_0);

            // Generate the skin model part
            ModelPart &skin_part = current_model.CreateModelPart("Boundaries");
            skin_part.CreateNewNode(662, 0.766593, 0.532174, 0.275516);
            skin_part.CreateNewNode(723, 0.793214, 0.506089, 0.308981);
            skin_part.CreateNewNode(737, 0.794158, 0.544627, 0.315665);
            skin_part.CreateNewNode(801, 0.81563, 0.518347, 0.349863);
            skin_part.CreateNewNode(814, 0.811818, 0.567485, 0.356072);
            skin_part.CreateNewNode(777, 0.809491, 0.469669, 0.339392);
            skin_part.CreateNewNode(710, 0.7901, 0.455512, 0.309309);
            skin_part.CreateNewNode(682, 0.768283, 0.578834, 0.289503);
            skin_part.CreateNewNode(741, 0.786372, 0.593624, 0.321883);
            skin_part.CreateNewNode(652, 0.766584, 0.482207, 0.273911);
            Properties::Pointer p_properties_1(new Properties(1));
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 477, {662,723,737}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 478, {737,723,801}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 479, {737,801,814}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 480, {801,723,777}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 510, {777,723,710}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 467, {682,662,737}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 484, {737,741,682}, p_properties_1);
            skin_part.CreateNewCondition("SurfaceCondition3D3N", 496, {723,652,710}, p_properties_1);

            // Call the intersections process
            FindIntersectedGeometricalObjectsProcess(volume_part, skin_part).Execute();

            // Check that there is no intersection
            KRATOS_CHECK(volume_part.GetElement(139).IsNot(SELECTED));
        }
    }
}  // namespace Kratos.
