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

#ifdef KRATOS_DEBUG
			KRATOS_DEBUG_CHECK(true);
			KRATOS_DEBUG_CHECK_IS_FALSE(false);

			KRATOS_DEBUG_CHECK_EQUAL(1.0, 1.0);
			KRATOS_DEBUG_CHECK_NOT_EQUAL(1.0, 2.0);

			KRATOS_DEBUG_CHECK_C_STRING_EQUAL("Test", "Test");
			KRATOS_DEBUG_CHECK_C_STRING_NOT_EQUAL("Test ", "Test");

			KRATOS_DEBUG_CHECK_LESS(1., 2.);
			KRATOS_DEBUG_CHECK_LESS_EQUAL(1., 1.);

			KRATOS_DEBUG_CHECK_GREATER(2., 1.);
			KRATOS_DEBUG_CHECK_GREATER_EQUAL(2., 2.);

			KRATOS_DEBUG_CHECK_STRING_CONTAIN_SUB_STRING(std::string("Test"), "es");
#else
			// all these checks should fail in debug mode but in release should be ignored
			KRATOS_DEBUG_CHECK(false); 
			KRATOS_DEBUG_CHECK_IS_FALSE(true);

			KRATOS_DEBUG_CHECK_EQUAL(1.0, 2.0);
			KRATOS_DEBUG_CHECK_NOT_EQUAL(1.0, 1.0);

			KRATOS_DEBUG_CHECK_C_STRING_EQUAL("T", "Test");
			KRATOS_DEBUG_CHECK_C_STRING_NOT_EQUAL("Test", "Test");

			KRATOS_DEBUG_CHECK_LESS(1., 1.);
			KRATOS_DEBUG_CHECK_LESS_EQUAL(1., 0.);

			KRATOS_DEBUG_CHECK_GREATER(1., 1.);
			KRATOS_DEBUG_CHECK_GREATER_EQUAL(0., 2.);

			KRATOS_DEBUG_CHECK_STRING_CONTAIN_SUB_STRING(std::string("Test"), "AnotherTEST");
#endif

		}
	}
}  // namespace Kratos.


