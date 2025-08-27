//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// Project includes
#include "testing/testing.h"
#include "future/processes/future_process.h"

namespace Kratos::Testing {

    /* This test will be added to a dedicated FutureCore Suit and will only be compiled with -DKRATOS_USE_FUTURE flag enabled.
     * For more info on future tests se also:
     * kratos/tests/cpp_tests/sources/test_namespaces.cpp
     */
    KRATOS_TEST_CASE_IN_SUITE(FutureProcess, KratosCoreFutureSuite)
    {
        Future::Process process;

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            process.Execute(),
            "I come from the future :)"
        );
    }
}