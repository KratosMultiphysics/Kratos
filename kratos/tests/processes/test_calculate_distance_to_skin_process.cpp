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
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_distance_to_skin_process.h"

namespace Kratos {
  namespace Testing {

	KRATOS_TEST_CASE_IN_SUITE(DistanceProcessQuadrilateral2D, KratosCoreFastSuite)
	{
		// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 10.00, 0.00, 0.00);

		Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

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
		CalculateDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

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

		KRATOS_CHECK_NEAR((surface_part.pGetNode(21))->FastGetSolutionStepValue(DISTANCE), -0.457143, 1e-6);
		KRATOS_CHECK_NEAR((surface_part.pGetNode(22))->FastGetSolutionStepValue(DISTANCE),  0.542857, 1e-6);
		KRATOS_CHECK_NEAR((surface_part.pGetNode(30))->FastGetSolutionStepValue(DISTANCE),  0.542857, 1e-6);
	}

	KRATOS_TEST_CASE_IN_SUITE(DistanceProcessSquareRing2D, KratosCoreFastSuite)
	{

		// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 10.00, 0.00, 0.00);

		Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

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
		CalculateDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

		KRATOS_CHECK_NEAR((surface_part.pGetNode(86))->FastGetSolutionStepValue(DISTANCE), 0.266667, 1e-6);
		KRATOS_CHECK_NEAR((surface_part.pGetNode(88))->FastGetSolutionStepValue(DISTANCE), -0.1, 1e-6);
		KRATOS_CHECK_NEAR((surface_part.pGetNode(112))->FastGetSolutionStepValue(DISTANCE), -0.566667, 1e-6);
		KRATOS_CHECK_NEAR((surface_part.pGetNode(138))->FastGetSolutionStepValue(DISTANCE),  0.733333, 1e-6);
	}

	KRATOS_TEST_CASE_IN_SUITE(HorizontalPlaneDistanceProcess, KratosCoreFastSuite)
	{

		// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point2 = Kratos::make_shared<Node<3>>(2, 10.00, 0.00, 0.00);
		Node<3>::Pointer p_point3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point4 = Kratos::make_shared<Node<3>>(4, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point5 = Kratos::make_shared<Node<3>>(5, 0.00, 0.00, 10.00);
		Node<3>::Pointer p_point6 = Kratos::make_shared<Node<3>>(6, 10.00, 0.00, 10.00);
		Node<3>::Pointer p_point7 = Kratos::make_shared<Node<3>>(7, 10.00, 10.00, 10.00);
		Node<3>::Pointer p_point8 = Kratos::make_shared<Node<3>>(8, 0.00, 10.00, 10.00);

		Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   1,
			"element_name":     "Element3D4N"
		})");

    Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(VELOCITY);
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
		StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

		// Generate the skin
		ModelPart &skin_part = current_model.CreateModelPart("Skin");
		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 0.0, 0.0, 2.0);
		skin_part.CreateNewNode(902, 10.0, 0.0, 2.0);
		skin_part.CreateNewNode(903, 10.0, 10.0, 2.0);
		skin_part.CreateNewNode(904, 0.0, 10.0, 2.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);

		// Compute distance
		CalculateDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

		for (auto& node : volume_part.Nodes())
			if (fabs(node.GetSolutionStepValue(DISTANCE)) < 1.00e16) { // There are no propagation in this version so I avoid numeric_limit::max() one
				auto distance = fabs(node.Z() - 2.00);
				KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), distance, 1e-6);
			}
	}

	KRATOS_TEST_CASE_IN_SUITE(HorizontalPlaneZeroDistanceProcess, KratosCoreFastSuite)
	{

		// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point2 = Kratos::make_shared<Node<3>>(2, 10.00, 0.00, 0.00);
		Node<3>::Pointer p_point3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point4 = Kratos::make_shared<Node<3>>(4, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point5 = Kratos::make_shared<Node<3>>(5, 0.00, 0.00, 10.00);
		Node<3>::Pointer p_point6 = Kratos::make_shared<Node<3>>(6, 10.00, 0.00, 10.00);
		Node<3>::Pointer p_point7 = Kratos::make_shared<Node<3>>(7, 10.00, 10.00, 10.00);
		Node<3>::Pointer p_point8 = Kratos::make_shared<Node<3>>(8, 0.00, 10.00, 10.00);

		Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   2,
			"element_name":     "Element3D4N"
		})");

    Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(VELOCITY);
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
		StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

		// Generate the skin
		ModelPart &skin_part = current_model.CreateModelPart("Skin");
		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 0.0, 0.0, 5.0);
		skin_part.CreateNewNode(902, 10.0, 0.0, 5.0);
		skin_part.CreateNewNode(903, 10.0, 10.0, 5.0);
		skin_part.CreateNewNode(904, 0.0, 10.0, 5.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,903,904 }, p_properties);

		// Compute distance
		CalculateDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

		for (auto& node : volume_part.Nodes())
			if (fabs(node.GetSolutionStepValue(DISTANCE)) < 1.00e16) { // There are no propagation in this version so I avoid numeric_limit::max() one
				auto distance = fabs(node.Z() - 5.00);
				KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), distance, 1e-6);
			}

	}

	KRATOS_TEST_CASE_IN_SUITE(TetrahedraInCubeDistanceProcess, KratosCoreFastSuite)
	{

		// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point2 = Kratos::make_shared<Node<3>>(2, 10.00, 0.00, 0.00);
		Node<3>::Pointer p_point3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point4 = Kratos::make_shared<Node<3>>(4, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point5 = Kratos::make_shared<Node<3>>(5, 0.00, 0.00, 10.00);
		Node<3>::Pointer p_point6 = Kratos::make_shared<Node<3>>(6, 10.00, 0.00, 10.00);
		Node<3>::Pointer p_point7 = Kratos::make_shared<Node<3>>(7, 10.00, 10.00, 10.00);
		Node<3>::Pointer p_point8 = Kratos::make_shared<Node<3>>(8, 0.00, 10.00, 10.00);

		Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

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
		CalculateDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

		KRATOS_CHECK_NEAR(volume_part.GetNode(135).GetSolutionStepValue(DISTANCE), 1.414213, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(136).GetSolutionStepValue(DISTANCE), 1.414213, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(137).GetSolutionStepValue(DISTANCE), 1.414213, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(256).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(257).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(258).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);

		//GidIO<> gid_io_fluid("C:/Temp/Tests/distance_test_fluid", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		//gid_io_fluid.InitializeMesh(0.00);
		//gid_io_fluid.WriteMesh(volume_part.GetMesh());
		//gid_io_fluid.FinalizeMesh();
		//gid_io_fluid.InitializeResults(0, volume_part.GetMesh());
		//gid_io_fluid.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
		//gid_io_fluid.FinalizeResults();

	}

	KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3IntersectionDistanceProcess, KratosCoreFastSuite)
	{

    Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.CreateNewNode(1, 1.00, 1.00, -10.00);
		volume_part.CreateNewNode(2, 1.00, 1.00, 10.00);
		volume_part.CreateNewNode(3, 10.00, 0.00, 0.00);
		volume_part.CreateNewNode(4, 0.00, 0.00, 0.00);

		Properties::Pointer p_properties(new Properties(0));
		volume_part.CreateNewElement("Element3D4N", 1, { 1,2,3,4 }, p_properties);

		// Generate the skin
		ModelPart &skin_part = current_model.CreateModelPart("Skin");
		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(11, 0.0, 0.0, 2.0);
		skin_part.CreateNewNode(12, 12.0, 0.0, 2.0);
		skin_part.CreateNewNode(13, 0.0, 12.0, 2.0);
		skin_part.CreateNewElement("Element3D3N", 1, { 11,12,13 }, p_properties);

		CalculateDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

		KRATOS_CHECK_NEAR(volume_part.GetNode(1).GetSolutionStepValue(DISTANCE), 12.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(2).GetSolutionStepValue(DISTANCE), 8.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(3).GetSolutionStepValue(DISTANCE), 2.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(4).GetSolutionStepValue(DISTANCE), 2.00, 1e-6);


	}

	KRATOS_TEST_CASE_IN_SUITE(Tetrahedra5IntersectionDistanceProcess, KratosCoreFastSuite)
	{

    Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.CreateNewNode(1, 2.50, 2.50, 0.00);
		volume_part.CreateNewNode(2, 2.50, 2.50, 2.50);
		volume_part.CreateNewNode(3, 2.50, 5.00, 2.50);
		volume_part.CreateNewNode(4, 5.00, 5.00, 2.50);

		Properties::Pointer p_properties(new Properties(0));
		volume_part.CreateNewElement("Element3D4N", 1, { 1,2,3,4 }, p_properties);

		// Generate the skin
		ModelPart &skin_part = current_model.CreateModelPart("Skin");
		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 6.0, 2.0, 2.0);
		skin_part.CreateNewNode(903, 4.0, 6.0, 2.0);
		skin_part.CreateNewNode(904, 4.0, 4.0, 7.0);

		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);

		CalculateDistanceToSkinProcess<3>(volume_part, skin_part).Execute();

		KRATOS_CHECK_NEAR(volume_part.GetNode(1).GetSolutionStepValue(DISTANCE), 2.0, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(2).GetSolutionStepValue(DISTANCE), -0.132068, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(3).GetSolutionStepValue(DISTANCE), 0.968496, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(4).GetSolutionStepValue(DISTANCE), 0.52827, 1e-6);
	}
}
}  // namespace Kratos.
