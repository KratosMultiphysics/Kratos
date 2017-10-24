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
//

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(FindIntersectedElementsProcess, KratosCoreFastSuite)
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

			ModelPart volume_part("Volume");
			ModelPart skin_part("Boundaries");
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
	}
}  // namespace Kratos.
