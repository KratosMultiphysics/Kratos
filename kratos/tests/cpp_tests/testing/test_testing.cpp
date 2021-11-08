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
#include <sstream>

// External includes

// Project includes
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(TestSuite, KratosCoreFastSuite) {
    Tester::CreateTestSuite("MyTestTestSuite");
    KRATOS_TESTING_ADD_TEST_TO_TEST_SUITE("TestTestSuite", "MyTestTestSuite");

    std::stringstream info;
    info << Tester::GetTestSuite("MyTestTestSuite");
    KRATOS_CHECK_NOT_EQUAL(
        info.str().find("MyTestTestSuite"), std::string::npos);
}
}
}  // namespace Kratos.
