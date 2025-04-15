//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Carlos A. Roig
//
//

#pragma once

// System includes

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "includes/expect.h"
#include "includes/data_communicator.h"

#include "testing/test_skipped_exception.h"

#include "tests/test_utilities/test_suite.h"            // Default Testing Suite. TODO: Remove once core is migrated to "core_test_suites".
#include "tests/test_utilities/test_environment.h"      // Environment used to initialize the tests.
#include "tests/test_utilities/test_event_listener.h"   // Custom Event Listener to control the output of the tests.

#define KRATOS_TEST_CASE(A) TEST_F(KratosCoreFastSuite, A)
#define KRATOS_TEST_CASE_IN_SUITE(A, B) TEST_F(B, A)
#define KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(A, B) TEST_F(B, A)

namespace Kratos::Testing {

/*
 * Initializes the parallel testing environment. This is usefull for other tests depending on a parallel environment.
*/
class GTestMain {
    public:
        static int InitializeTesting(int argc, char* argv[]) {
            // Initialize the tests
            ::testing::InitGoogleTest(&argc, argv);

            // Remove the default listener
            testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
            auto default_printer = listeners.Release(listeners.default_result_printer());

            // Create a configurable listener
            Kratos::Testing::ConfigurableEventListener *listener = new Kratos::Testing::ConfigurableEventListener(default_printer);

            // Add the common environment to the test
            ::testing::AddGlobalTestEnvironment(new Kratos::Testing::KratosTestEnv);

            // Add our listener
            listeners.Append(listener);

            // Run the tests
            return RUN_ALL_TESTS();
        }
};

} // namespace Kratos::Testing
