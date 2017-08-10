//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/check_utilities.h"
#include "includes/model_part.h"

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(CheckUtilities, KratosCoreFastSuite)
        {
			ModelPart model_part("TestModelPart");

            model_part.AddNodalSolutionStepVariable(VELOCITY);
            model_part.AddNodalSolutionStepVariable(PRESSURE);

            model_part.SetBufferSize(1);

            model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

            for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
            {
                it_node->AddDof(VELOCITY_Y,REACTION_Y);
                it_node->AddDof(PRESSURE);
            }

            Node<3> & r_node = *(model_part.NodesBegin());

            // These functions throw an error if the check fails

            // Expected passes: test OK if no error is trown
            CheckUtilities::CheckVariableKey(VELOCITY);
            CheckUtilities::CheckVariableInNodalData(VELOCITY, r_node);
            CheckUtilities::CheckVariableInNodalData(PRESSURE, r_node);
            CheckUtilities::CheckVariableInNodalData(VELOCITY_X, r_node);
            CheckUtilities::CheckDofInNode(VELOCITY_Y, r_node);
            CheckUtilities::CheckDofInNode(PRESSURE, r_node);

            // Expected fails: test failed unless an error is thrown
            KRATOS_CHECK_ERROR_IS_THROWN( CheckUtilities::CheckVariableInNodalData(BODY_FORCE, r_node), "Missing BODY_FORCE variable in solution step data for node 1." );
            KRATOS_CHECK_ERROR_IS_THROWN( CheckUtilities::CheckVariableInNodalData(BODY_FORCE_X, r_node), "Missing BODY_FORCE_X variable in solution step data for node 1." );
            KRATOS_CHECK_ERROR_IS_THROWN( CheckUtilities::CheckVariableInNodalData(TEMPERATURE, r_node), "Missing TEMPERATURE variable in solution step data for node 1." );
            KRATOS_CHECK_ERROR_IS_THROWN( CheckUtilities::CheckDofInNode(TEMPERATURE, r_node), "Missing Degree of Freedom for TEMPERATURE in node 1." );
            KRATOS_CHECK_ERROR_IS_THROWN( CheckUtilities::CheckDofInNode(VELOCITY_X, r_node), "Missing Degree of Freedom for VELOCITY_X in node 1." );
        }
    }
}


