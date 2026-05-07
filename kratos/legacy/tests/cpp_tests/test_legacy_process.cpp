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
#include "legacy/processes/legacy_process.h"

namespace Kratos::Testing {

    /* This test will be added to a dedicated LegacyCore Suit and will only be compiled with -DKRATOS_USE_LEGACY flag enabled.
     * For more info on legacy tests se also:
     * kratos/tests/cpp_tests/sources/test_namespaces.cpp
     */
    KRATOS_TEST_CASE_IN_SUITE(LegacyProcess, KratosCoreLegacySuite)
    {
        Legacy::Process process;

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            process.Execute(),
            "I... I don't feel well... My time has come, soon I will be unmade. Thank you for all the segmentation faults we lived together. :_)"
        );
    }
}