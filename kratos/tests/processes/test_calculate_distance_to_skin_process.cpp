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
#include "includes/checks.h"
#include "includes/gid_io.h"
// #include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/calculate_signed_distance_to_3d_skin_process.h" // TODO: Change the tested process as soon as the new distance process is available
#include "processes/structured_mesh_generator_process.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(CalculateDistanceToSkinProcess, KratosCoreFastSuite)
		{

			// Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
			Point<3>::Pointer p_point1(new Point<3>(0.00, 0.00, 0.00));
			Point<3>::Pointer p_point2(new Point<3>(10.00, 0.00, 0.00));
			Point<3>::Pointer p_point3(new Point<3>(10.00, 10.00, 0.00));
			Point<3>::Pointer p_point4(new Point<3>(0.00, 10.00, 0.00));
			Point<3>::Pointer p_point5(new Point<3>(0.00, 0.00, 10.00));
			Point<3>::Pointer p_point6(new Point<3>(10.00, 0.00, 10.00));
			Point<3>::Pointer p_point7(new Point<3>(10.00, 10.00, 10.00));
			Point<3>::Pointer p_point8(new Point<3>(0.00, 10.00, 10.00));

			Hexahedra3D8<Point<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

			Parameters mesher_parameters(R"(
            {
                "number_of_divisions": 	20,
                "element_name": 		"Element3D4N"
            })");

			ModelPart volume_part("Volume");
			volume_part.AddNodalSolutionStepVariable(VELOCITY);
			volume_part.AddNodalSolutionStepVariable(DISTANCE);
			volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
			StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

			// Generate the skin
			ModelPart skin_part("Skin");
			skin_part.AddNodalSolutionStepVariable(VELOCITY);
			skin_part.CreateNewNode(1, 2.0, 2.0, 2.0);
			skin_part.CreateNewNode(2, 6.0, 2.0, 2.0);
			skin_part.CreateNewNode(3, 4.0, 6.0, 2.0);
			// skin_part.CreateNewNode(4, 4.0, 6.0, 7.0);
			skin_part.CreateNewNode(4, 4.0, 4.0, 7.0);
			Properties::Pointer p_properties(new Properties(0));
			skin_part.CreateNewElement("Element3D3N", 1, { 1,2,3 }, p_properties);
			skin_part.CreateNewElement("Element3D3N", 2, { 1,4,3 }, p_properties);
			skin_part.CreateNewElement("Element3D3N", 3, { 2,3,4 }, p_properties);
			skin_part.CreateNewElement("Element3D3N", 4, { 1,2,4 }, p_properties);

			// Compute distance
			// TODO: Change the tested process as soon as the new distance process is available
			CalculateSignedDistanceTo3DSkinProcess(skin_part, volume_part).Execute();

			// GidIO<> gid_io_fluid("/home/rzorrilla/Kratos/tests/distance_test_fluid", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
			// gid_io_fluid.InitializeMesh(0.00);
			// gid_io_fluid.WriteMesh(volume_part.GetMesh());
			// gid_io_fluid.FinalizeMesh();
			// gid_io_fluid.InitializeResults(0, volume_part.GetMesh());
			// gid_io_fluid.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
			// gid_io_fluid.FinalizeResults();
			//
			// GidIO<> gid_io_skin("/home/rzorrilla/Kratos/tests/distance_test_skin", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
			// gid_io_skin.InitializeMesh(0.00);
			// gid_io_skin.WriteMesh(skin_part.GetMesh());
			// gid_io_skin.FinalizeMesh();

		}
	}
}  // namespace Kratos.
