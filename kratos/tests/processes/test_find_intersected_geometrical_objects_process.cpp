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
                "number_of_divisions":2,
                "element_name": "Element3D4N"
            }  )");

			std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

			ModelPart volume_part("Volume");
			ModelPart skin_part("Boundaries");
			skin_part.CreateNewNode(1, 1., 1., 1.);
			skin_part.CreateNewNode(2, 1., 2., 2.);
			skin_part.CreateNewNode(3, 4., 1., 1.);
			Properties::Pointer p_properties(new Properties(0));
			skin_part.CreateNewElement("Element3D3N", 1, { 1,2,3 }, p_properties);
			StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();
			FindIntersectedGeometricalObjectsProcess(volume_part, skin_part).Execute();
			std::cout << std::endl;
			for (auto& element : volume_part.Elements()) {
				std::cout << element.Id() << " : " << element.Is(SELECTED) << std::endl;
			}
		}
	}
}  // namespace Kratos.
