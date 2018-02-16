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
#include "includes/prime_numbers.h"
#include "includes/model_part.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(NodeAssignOperator, KratosCoreFastSuite)
		{
			ModelPart model_part("test");
			model_part.AddNodalSolutionStepVariable(DISTANCE);
			model_part.AddNodalSolutionStepVariable(VELOCITY);
            
			auto p_node = model_part.CreateNewNode(1, 1., 0, 0);

            p_node->FastGetSolutionStepValue(DISTANCE) = 12.1;
            p_node->FastGetSolutionStepValue(VELOCITY_X) = 32.4;
            p_node->Set(ACTIVE, true);

            Node<3> copy_of_node(2,1,0,0);
            copy_of_node = *p_node;

            KRATOS_CHECK_EQUAL(copy_of_node.Id(), 1);
            KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(DISTANCE), 12.1);
            KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_X), 32.4);
            KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_Y), 0.00);
            KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_Z), 0.00);
            KRATOS_CHECK(copy_of_node.Is(ACTIVE));
		}
	}
}  // namespace Kratos.
