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

namespace Kratos {
  namespace Testing {


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorTetrahedraSkinProcess, KratosCoreFastSuite)
	{
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
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

	Point dummy(0.00,0.00,0.00);
	for(auto& node : volume_part.Nodes()){
		if(tetrahedra.IsInside(node.Coordinates(),dummy)){
			KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
		}
	}
		// Note that we cannot check the outside because on the interface is not well defined
		KRATOS_CHECK_NEAR(volume_part.GetNode(135).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(136).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
		KRATOS_CHECK_NEAR(volume_part.GetNode(137).GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
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


}
}  // namespace Kratos.
