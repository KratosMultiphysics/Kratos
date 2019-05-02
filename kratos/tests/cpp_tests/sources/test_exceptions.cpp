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

		KRATOS_TEST_CASE_IN_SUITE(ExceptionDefaultConstruction, KratosCoreFastSuite)
		{
			try {
				throw Exception();
			}
			catch (Exception& e) {
				KRATOS_CHECK_C_STRING_EQUAL(e.what(), "Unknown Error\nin Unknown Location");
				KRATOS_CHECK_EQUAL(e.where().CleanFileName(), "Unknown File");
				KRATOS_CHECK_EQUAL(e.where().CleanFunctionName(), "Unknown Location");
				KRATOS_CHECK_EQUAL(e.where().GetLineNumber(), 0);
			}
		}


	}
}  // namespace Kratos.
