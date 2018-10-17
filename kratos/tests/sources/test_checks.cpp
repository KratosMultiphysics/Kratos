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

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(Checks, KratosCoreFastSuite)
		{
			KRATOS_CHECK(true);
			KRATOS_CHECK_IS_FALSE(false);

			KRATOS_CHECK_EQUAL(1.0, 1.0);
			KRATOS_CHECK_NOT_EQUAL(1.0, 2.0);

			KRATOS_CHECK_C_STRING_EQUAL("Test", "Test");
			KRATOS_CHECK_C_STRING_NOT_EQUAL("Test ", "Test");

			KRATOS_CHECK_LESS(1., 2.);
			KRATOS_CHECK_LESS_EQUAL(1., 1.);

			KRATOS_CHECK_GREATER(2., 1.);
			KRATOS_CHECK_GREATER_EQUAL(2., 2.);

			KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(std::string("Test"), "es");
		}

		KRATOS_TEST_CASE_IN_SUITE(VariableChecks, KratosCoreFastSuite)
        {
            Model current_model;
            
            ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

            model_part.AddNodalSolutionStepVariable(VELOCITY);
            model_part.AddNodalSolutionStepVariable(PRESSURE);
            model_part.AddNodalSolutionStepVariable(REACTION);

            model_part.SetBufferSize(1);

            model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

            for (auto it_node = model_part.NodesBegin();
                 it_node != model_part.NodesEnd(); ++it_node)
            {
                it_node->AddDof(VELOCITY_Y, REACTION_Y);
                it_node->AddDof(PRESSURE);
            }

            Node<3>& r_node = *(model_part.NodesBegin());

            // These functions throw an error if the check fails
            // Expected passes: test OK if no error is thrown
            KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);

            // Expected fails: test failed unless an error is thrown
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node),
                "Missing BODY_FORCE variable in solution step data for node "
                "1.");
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE_X, r_node),
                "Missing BODY_FORCE_X variable in solution step data for node "
                "1.");
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE, r_node),
                "Missing TEMPERATURE variable in solution step data for node "
                "1.");
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_DOF_IN_NODE(TEMPERATURE, r_node),
                "Missing Degree of Freedom for TEMPERATURE in node 1.");
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node),
                "Missing Degree of Freedom for VELOCITY_X in node 1.");
        }
    }
}  // namespace Kratos.
