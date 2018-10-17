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
#include "geometries/hexahedra_3d_8.h"
#include "includes/gid_io.h"


namespace Kratos {
	namespace Testing {


		KRATOS_TEST_CASE_IN_SUITE(TetrahedraMeshEdgeSwappingProcess, KratosCoreFastSuite)
		{
			Model current_model;
			
			Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
			Node<3>::Pointer p_point2(new Node<3>(2, 10.00, 0.00, 0.00));
			Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
			Node<3>::Pointer p_point4(new Node<3>(4, 0.00, 10.00, 0.00));
			Node<3>::Pointer p_point5(new Node<3>(5, 0.00, 0.00, 10.00));
			Node<3>::Pointer p_point6(new Node<3>(6, 10.00, 0.00, 10.00));
			Node<3>::Pointer p_point7(new Node<3>(7, 10.00, 10.00, 10.00));
			Node<3>::Pointer p_point8(new Node<3>(8, 0.00, 10.00, 10.00));

			Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

			ModelPart& model_part = current_model.CreateModelPart("Test");

			Parameters mesher_parameters(R"(
            {
                "number_of_divisions":1,
                "element_name": "Element3D4N"
            }  )");

			//std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

			StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

			GidIO<> gid_io("/home/pooyan/kratos/tests/edge_swapping", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
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

	}
}  // namespace Kratos.
