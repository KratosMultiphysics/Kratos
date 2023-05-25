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
#include <vector>

// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
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

            Node& r_node = *(model_part.NodesBegin());

            // These functions throw an error if the check fails
            // Expected passes: test OK if no error is thrown
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

        KRATOS_TEST_CASE_IN_SUITE(VectorChecks, KratosCoreFastSuite)
        {
            std::vector<double> v1{1.0, 2.0};
            std::vector<double> v2{1.0, 1.99};
            std::vector<double> v3{1.0, 2.0, 3.0};

            KRATOS_CHECK_VECTOR_EQUAL(v1, v1);
            KRATOS_CHECK_VECTOR_NEAR(v1, v2, 0.1);

            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_EQUAL(v1, v3),
                "Check failed because vector arguments do not have the same size"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_NEAR(v1, v3, 0.1),
                "Check failed because vector arguments do not have the same size"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_EQUAL(v1,v2),
                "Mismatch found in component 1:"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_NEAR(v1,v2, 0.001),
                "Mismatch found in component 1:"
            );
        }

        KRATOS_TEST_CASE_IN_SUITE(MatrixChecks, KratosCoreFastSuite)
        {
            Matrix m1 = IdentityMatrix(2);
            Matrix m2 = IdentityMatrix(2);
            m2(0,1) = 0.01;
            Matrix m3 = IdentityMatrix(3);

            KRATOS_CHECK_MATRIX_EQUAL(m1, m1);
            KRATOS_CHECK_MATRIX_NEAR(m1, m2, 0.1);

            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_EQUAL(m1, m3),
                "Check failed because matrix arguments do not have the same dimensions"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_NEAR(m1, m3, 0.1),
                "Check failed because matrix arguments do not have the same dimensions"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_EQUAL(m1,m2),
                "Mismatch found in component (0,1):"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_NEAR(m1,m2, 0.001),
                "Mismatch found in component (0,1):"
            );
        }

        KRATOS_TEST_CASE_IN_SUITE(NaNValuesCheck, KratosCoreFastSuite)
		{
            double a = std::numeric_limits<double>::quiet_NaN();
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_NEAR(1.0, a, 1e-7),
                "Check failed because 1.0 = 1 is not near to a = nan within the tolerance 1e-07"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_NEAR(a, 1.0, 1e-7),
                "Check failed because a = nan is not near to 1.0 = 1 within the tolerance 1e-07"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_LESS_EQUAL(a,1.0),
                "Check failed because a is greater than 1.0"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_LESS_EQUAL(1.0,a),
                "Check failed because 1.0 is greater than a"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_GREATER_EQUAL(a,1.0),
                "Check failed because a is less than 1.0"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_GREATER_EQUAL(1.0,a),
                "Check failed because 1.0 is less than a"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_RELATIVE_NEAR(1.0,a,1e-3),
                "Check failed because 1.0 = 1 is not near to a = nan within the relative tolerance 0.001"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_RELATIVE_NEAR(a,1.0,1e-3),
                "Check failed because a = nan is not near to 1.0 = 1 within the relative tolerance 0.001"
            );

            std::vector<double> v1{1.0, 2.0};
            std::vector<double> v2{1.0, a};
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_NEAR(v1,v2,1e-3),
                "Check failed because vector v1 with values\n[1, 2]\n"
                "Is not near vector v2 with values\n[1, nan]"
                "\nMismatch found in component 1:\n2 not near nan within tolerance 0.001."
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_NEAR(v2,v1,1e-3),
                "Check failed because vector v2 with values\n[1, nan]\n"
                "Is not near vector v1 with values\n[1, 2]"
                "\nMismatch found in component 1:\nnan not near 2 within tolerance 0.001."
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_RELATIVE_NEAR(v1,v2,1e-3),
                "Check failed because vector v1 with values\n[1, 2]\n"
                "Is not near vector v2 with values\n[1, nan]"
                "\nMismatch found in component 1:\n2 not near nan within tolerance 0.001."
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_VECTOR_RELATIVE_NEAR(v2,v1,1e-3),
                "Check failed because vector v2 with values\n[1, nan]\n"
                "Is not near vector v1 with values\n[1, 2]"
                "\nMismatch found in component 1:\nnan not near 2 within relative tolerance 0.001."
            );

            Matrix m1 = IdentityMatrix(2);
            Matrix m2 = IdentityMatrix(2);
            m2(0,1) = a;
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_NEAR(m1,m2,1e-3),
                "Check failed because matrix m1 with values\n[2,2]((1,0),(0,1))\n"
                "Is not near matrix m2 with values\n[2,2]((1,nan),(0,1))"
                "\nMismatch found in component (0,1): \n0 not near nan within tolerance 0.001.\n"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_NEAR(m2,m1,1e-3),
                "Check failed because matrix m2 with values\n[2,2]((1,nan),(0,1))\n"
                "Is not near matrix m1 with values\n[2,2]((1,0),(0,1))\n"
                "Mismatch found in component (0,1): \nnan not near 0 within tolerance 0.001.\n"
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_RELATIVE_NEAR(m1,m2,1e-3),
                "Check failed because matrix m1 with values\n[2,2]((1,0),(0,1))\n"
                "Is not near matrix m2 with values\n[2,2]((1,nan),(0,1))\n"
                "Mismatch found in component (0,1): \n0 not near nan within tolerance 0.001."
            );
            KRATOS_CHECK_EXCEPTION_IS_THROWN(
                KRATOS_CHECK_MATRIX_RELATIVE_NEAR(m2,m1,1e-3),
                "Check failed because matrix m2 with values\n[2,2]((1,nan),(0,1))\n"
                "Is not near matrix m1 with values\n[2,2]((1,0),(0,1))\n"
                "Mismatch found in component (0,1): \nnan not near 0 within tolerance 0.001."
            );
		}
    }
}  // namespace Kratos.
