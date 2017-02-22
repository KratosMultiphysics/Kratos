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
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/tetrahedra_mesh_edge_swapping_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/gid_io.h"


namespace Kratos {
	namespace Testing {


		KRATOS_TEST_CASE_IN_SUITE(TetrahedraMeshEdgeSwappingProcess, KratosCoreFastSuite)
		{
			Point<3>::Pointer p_point1(new Point<3>(0.00, 0.00, 0.00));
			Point<3>::Pointer p_point2(new Point<3>(10.00, 0.00, 0.00));
			Point<3>::Pointer p_point3(new Point<3>(10.00, 10.00, 0.00));
			Point<3>::Pointer p_point4(new Point<3>(0.00, 10.00, 0.00));
			Point<3>::Pointer p_point5(new Point<3>(0.00, 0.00, 10.00));
			Point<3>::Pointer p_point6(new Point<3>(10.00, 0.00, 10.00));
			Point<3>::Pointer p_point7(new Point<3>(10.00, 10.00, 10.00));
			Point<3>::Pointer p_point8(new Point<3>(0.00, 10.00, 10.00));

			Hexahedra3D8<Point<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

			ModelPart model_part("Test");

			Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element3D4N"
            }  )");

			//std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

			StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

			GidIO<> gid_io("c:/temp/coarsening/edge_swapping_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
			gid_io.InitializeMesh(0.00);
			gid_io.WriteMesh(model_part.GetMesh());
			gid_io.FinalizeMesh();

			ModelPart& r_skin_model_part = model_part.GetSubModelPart("Skin");

			for (auto i_node = r_skin_model_part.NodesBegin(); i_node != r_skin_model_part.NodesEnd(); i_node++)
				i_node->Set(BOUNDARY);

			FindNodalNeighboursProcess(model_part).Execute();

			TetrahedraMeshEdgeSwappingProcess(model_part).Execute();

			// gid_io.InitializeMesh(1.00);
			// gid_io.WriteMesh(model_part.GetMesh());
			// gid_io.FinalizeMesh();



		}
		
		KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3to2EdgeSwappingProcess, KratosCoreFastSuite)
		{
			ModelPart model_part("Test");

			model_part.CreateNewNode(1, 0.00, 0.00, 0.00);
			model_part.CreateNewNode(2, 10.00, 0.00, 0.00);
			model_part.CreateNewNode(3, 10.00, 10.00, 0.00);
			model_part.CreateNewNode(4, 1.00, 1.00, -10.00);
			model_part.CreateNewNode(5, 1.00, 1.00, 10.00);

			Properties::Pointer p_properties(new Properties(0));
			model_part.CreateNewElement("Element3D4N", 1, { 4,5,3,1 }, p_properties);
			model_part.CreateNewElement("Element3D4N", 2, { 4,5,2,3 }, p_properties);
			model_part.CreateNewElement("Element3D4N", 3, { 5,4,2,1 }, p_properties);


			FindNodalNeighboursProcess(model_part).Execute();

			GidIO<> gid_io("c:/temp/coarsening/edge_swapping_3to2_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
			gid_io.InitializeMesh(0.00);
			gid_io.WriteMesh(model_part.GetMesh());
			gid_io.FinalizeMesh();
			TetrahedraMeshEdgeSwappingProcess(model_part).Execute();

			KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), 2);
			KRATOS_CHECK_GREATER(model_part.GetElement(1).GetGeometry().Volume(), 166.);
			KRATOS_CHECK_GREATER(model_part.GetElement(2).GetGeometry().Volume(), 166.);

			 gid_io.InitializeMesh(1.00);
			 gid_io.WriteMesh(model_part.GetMesh());
			 gid_io.FinalizeMesh();



		}

	}
}  // namespace Kratos.
