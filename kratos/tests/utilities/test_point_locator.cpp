//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/point_locator.h"

namespace Kratos {
    namespace Testing {


        KRATOS_TEST_CASE_IN_SUITE(PointLocatorTriangle, KratosCoreFastSuite)
        {
            // Generate a model part with the previous
			ModelPart model_part("Triangle");
            model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

			// Fill the model part geometry data
			model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			Properties::Pointer p_properties(new Properties(0));
            const int elem_id = 43;
			model_part.CreateNewElement("Element2D3N", elem_id, {1, 2, 3}, p_properties);

            PointLocator point_locator(model_part);

            Point the_point(0.1, 0.25, 0.0);

            int found_id = -1;

            Vector shape_function_values;

            const bool search_successful = point_locator.FindElement(the_point, found_id, shape_function_values);

            KRATOS_CHECK(search_successful);
            KRATOS_CHECK_EQUAL(found_id, elem_id);
            KRATOS_CHECK_EQUAL(shape_function_values.size(), 3);
            KRATOS_CHECK_NEAR(shape_function_values[0], 0.65, 1e-12);
            KRATOS_CHECK_NEAR(shape_function_values[1], 0.1,  1e-12);
            KRATOS_CHECK_NEAR(shape_function_values[2], 0.25, 1e-12);
            // std::cout << shape_function_values[0] << " "<<shape_function_values[1] << " "<<shape_function_values[2] << std::endl;
        }

        KRATOS_TEST_CASE_IN_SUITE(PointLocatorQuadrilateral, KratosCoreFastSuite)
        {
        }

        KRATOS_TEST_CASE_IN_SUITE(PointLocatorTetrahedra, KratosCoreFastSuite)
        {
        }

        KRATOS_TEST_CASE_IN_SUITE(PointLocatorHexahedra, KratosCoreFastSuite)
        {
        }

    } // namespace Testing
}  // namespace Kratos.
