//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
    #include "future/processes/future_process.h"
#endif
#ifdef KRATOS_USE_LEGACY
    #include "legacy/processes/legacy_process.h"
#endif

namespace Kratos::Testing
{

/* Tests in this file will be always compiled and will change behaviour based on the value of -DKRATOS_USE_FUTURE and -DKRATOS_USE_LEGACY
 * For tests for specific codes that only exist inside one of the namespaces see:
 * - kratos/future/tests/cpp_tests/test_future_process.cpp
 * - kratos/legacy/tests/cpp_tests/test_legacy_process.cpp
 */

/* Test if the class exists in the namespace if the namespace has been compiled */
KRATOS_TEST_CASE_IN_SUITE(FutureProcess, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_FUTURE
        Future::Process process;

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            process.Execute(),
            "I come from the future :)"
        );
    #else
        KRATOS_SUCCEED();
    #endif
}

/* Test if the class exists in the registry if the namespace has been compiled */
KRATOS_TEST_CASE_IN_SUITE(FutureProcessFromRegistry, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_FUTURE
        KRATOS_EXPECT_TRUE(Registry::HasItem("Processes.KratosMultiphysics.Future.Process"))

        if (Registry::HasItem("Processes.KratosMultiphysics.Future")) {
            auto process = Registry::GetValue<Future::Process>("Processes.KratosMultiphysics.Future.Process.Prototype");

            KRATOS_EXPECT_EXCEPTION_IS_THROWN(
                process.Execute(),
                "I come from the future :)"
            );
        }
    #else
        KRATOS_EXPECT_FALSE(Registry::HasItem("Processes.KratosMultiphysics.Future.Process"))
    #endif
}

/* Test that the class in the registry is not the same class as its legacy / production counterpart */
KRATOS_TEST_CASE_IN_SUITE(FutureProcessFromRegistryIsNotProcess, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_FUTURE
        KRATOS_EXPECT_TRUE(Registry::HasItem("Processes.KratosMultiphysics.Future.Process"))

        if (Registry::HasItem("Processes.KratosMultiphysics.Future")) {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(
                auto process = Registry::GetValue<Process>("Processes.KratosMultiphysics.Future.Process.Prototype");,
                "Error: bad any_cast"
            );
        }
    #else
        KRATOS_EXPECT_FALSE(Registry::HasItem("Processes.KratosMultiphysics.Future.Process"))
    #endif
}

/* Test if the class exists in the namespace if the namespace has been compiled */
KRATOS_TEST_CASE_IN_SUITE(LegacyProcess, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_LEGACY
        Legacy::Process process;

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            process.Execute(),
            "I... I don't feel well... My time has come, soon I will be unmade. Thank you for all the segmentation faults we lived together. :_)"
        );
    #else
        KRATOS_SUCCEED();
    #endif
}

/* Test if the class exists in the registry if the namespace has been compiled */
KRATOS_TEST_CASE_IN_SUITE(LegacyProcessFromRegistry, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_LEGACY
        KRATOS_EXPECT_TRUE(Registry::HasItem("Processes.KratosMultiphysics.Legacy.Process"))

        if (Registry::HasItem("Processes.KratosMultiphysics.Legacy.Process")) {
            auto process = Registry::GetValue<Legacy::Process>("Processes.KratosMultiphysics.Legacy.Process.Prototype");

            KRATOS_EXPECT_EXCEPTION_IS_THROWN(
                process.Execute(),
                "I... I don't feel well... My time has come, soon I will be unmade. Thank you for all the segmentation faults we lived together. :_)"
            );
        }
    #else
        KRATOS_EXPECT_FALSE(Registry::HasItem("Processes.KratosMultiphysics.Legacy.Process"))
    #endif
}

/* Test that the class in the registry is not the same class as its future / production counterpart */
KRATOS_TEST_CASE_IN_SUITE(LegacyProcessFromRegistryIsNotProcess, KratosCoreFastSuite)
{
    #ifdef KRATOS_USE_LEGACY
        KRATOS_EXPECT_TRUE(Registry::HasItem("Processes.KratosMultiphysics.Legacy.Process"))

        if (Registry::HasItem("Processes.KratosMultiphysics.Legacy")) {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(
                auto process = Registry::GetValue<Process>("Processes.KratosMultiphysics.Legacy.Process.Prototype");,
                "Error: bad any_cast"
            );
        }
    #else
        KRATOS_EXPECT_FALSE(Registry::HasItem("Processes.KratosMultiphysics.Legacy.Process"))
    #endif
}

} // Kratos::Testing namespace
