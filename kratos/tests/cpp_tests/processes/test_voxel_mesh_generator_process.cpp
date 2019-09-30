//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//					 Ruben Zorrilla
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
// #include "includes/gid_io.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "processes/voxel_mesh_generator_process.h"
#include "processes/voxel_mesh_coloring_process.h"

namespace Kratos {
  namespace Testing {


	KRATOS_TEST_CASE_IN_SUITE(XCartesianRayplaneIntersection, KratosCoreFastSuite)
	{
		Point::Pointer p_point_1=Kratos::make_shared<Point>(.4, 0.00, 0.00);
		Point::Pointer p_point_2=Kratos::make_shared<Point>(.4, 1.00, 0.00);
		Point::Pointer p_point_3=Kratos::make_shared<Point>(.4, 1.00, 1.00);
		Point::Pointer p_point_4=Kratos::make_shared<Point>(.4, 0.00, 1.00);

		Triangle3D3<Point> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Point> triangle_2(p_point_1, p_point_3, p_point_4);

		Kratos::Internals::CartesianRay ray_1(0, Point(0.00, 0.50, 0.50), Point(1.00, 0.50, 0.50));
		Kratos::Internals::CartesianRay ray_2(0, Point(0.00, 0.00, 0.00), Point(1.00, 0.00, 0.00));
		Kratos::Internals::CartesianRay ray_3(0, Point(0.00, 0.20, 0.70), Point(1.00, 0.20, 0.70));

		ray_1.AddIntersection(triangle_1, 1e-9);
		ray_2.AddIntersection(triangle_1, 1e-9);
		ray_3.AddIntersection(triangle_1, 1e-9);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints().size(), 1);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersectionPoints().size(), 0);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints()[0], .4);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints()[0], .4);

		ray_1.AddIntersection(triangle_2, 1e-9);
		ray_2.AddIntersection(triangle_2, 1e-9);
		ray_3.AddIntersection(triangle_2, 1e-9);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints().size(), 2);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints().size(), 2);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersectionPoints().size(), 1);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints()[1], .4);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints()[1], .4);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersectionPoints()[0], .4);

		ray_1.CollapseIntersectionPoints(1e-9);
		ray_2.CollapseIntersectionPoints(1e-9);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints().size(), 1);
		
		KRATOS_CHECK_EQUAL(ray_1.GetIntersectionPoints()[0], .4);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersectionPoints()[0], .4);
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessNodesPositions, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [5,5,5],
			"element_name":     "Element3D4N"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{1.00, 2.00, 3.00}, Point{11.00, 12.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k <= 5; k++) {
			for (std::size_t j = 0; j <= 5; j++) {
				for (std::size_t i = 0; i <= 5; i++) {
					auto& node = *i_node++;
					double x = 2.00*i + 1.00;
					double y = 2.00*j + 2.00;
					double z = 2.00*k + 3.00;
                	KRATOS_CHECK_NEAR(node.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z, 1e-6);
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCentersPositions, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{1.00, 2.00, 3.00}, Point{11.00, 12.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k < 10; k++) {
			for (std::size_t j = 0; j < 10; j++) {
				for (std::size_t i = 0; i < 10; i++) {
					auto& node = *i_node++;
					double x = i+1.50;
					double y = j+2.50;
					double z = k+3.50;
                	KRATOS_CHECK_NEAR(node.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z, 1e-6);
				}
            }
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessTetrahedraElementColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(VELOCITY);
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

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

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();
		// Compute distance
		// VoxelMeshColoringProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

        Point dummy(0.00,0.00,0.00);
        for(auto& element : volume_part.Elements()){
            if(tetrahedra.IsInside(element.GetGeometry().Center(),dummy)){
                KRATOS_CHECK_NEAR(element.GetValue(DISTANCE), -1.00, 1e-6);
            }
		}

		//GidIO<> gid_io_fluid("C:/Temp/Tests/distance_test_fluid", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		//gid_io_fluid.InitializeMesh(0.00);
		//gid_io_fluid.WriteMesh(volume_part.GetMesh());
		//gid_io_fluid.FinalizeMesh();
		//gid_io_fluid.InitializeResults(0, volume_part.GetMesh());
		//gid_io_fluid.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
		//gid_io_fluid.FinalizeResults();

	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessTetrahedraElementNodesColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "nodes"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(VELOCITY);
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

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

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();
		// Compute distance
		// VoxelMeshColoringProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

        Point dummy(0.00,0.00,0.00);
        for(auto& node : volume_part.Nodes()){
            if(!tetrahedra.IsInside(node,dummy)){
				// KRATOS_WATCH(node.GetSolutionStepValue(DISTANCE));
                KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
            }
		}
	}


}
}  // namespace Kratos.
