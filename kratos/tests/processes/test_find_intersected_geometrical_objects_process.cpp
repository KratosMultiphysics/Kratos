//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
// #include "includes/gid_io.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcess2D, KratosCoreFastSuite)
		{     
			Node<3>::Pointer p_point1(new Node<3>(1, 0.0, 0.0, 0.0));
			Node<3>::Pointer p_point2(new Node<3>(2, 0.0, 1.0, 0.0));
			Node<3>::Pointer p_point3(new Node<3>(3, 1.0, 1.0, 0.0));
			Node<3>::Pointer p_point4(new Node<3>(4, 1.0, 0.0, 0.0));
      
			Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);
      
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

		KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcessNoIntersection2D, KratosCoreFastSuite)
		{
			Node<3>::Pointer p_point1(new Node<3>(1, 0.0, 0.0, 0.0));
			Node<3>::Pointer p_point2(new Node<3>(2, 0.0, 1.0, 0.0));
			Node<3>::Pointer p_point3(new Node<3>(3, 1.0, 1.0, 0.0));
			Node<3>::Pointer p_point4(new Node<3>(4, 1.0, 0.0, 0.0));

			Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

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

		KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcess3D, KratosCoreFastSuite)
		{
			Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
			Node<3>::Pointer p_point2(new Node<3>(2, 10.00, 0.00, 0.00));
			Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
			Node<3>::Pointer p_point4(new Node<3>(4, 0.00, 10.00, 0.00));
			Node<3>::Pointer p_point5(new Node<3>(5, 0.00, 0.00, 10.00));
			Node<3>::Pointer p_point6(new Node<3>(6, 10.00, 0.00, 10.00));
			Node<3>::Pointer p_point7(new Node<3>(7, 10.00, 10.00, 10.00));
			Node<3>::Pointer p_point8(new Node<3>(8, 0.00, 10.00, 10.00));

			Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

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
	}
}  // namespace Kratos.
