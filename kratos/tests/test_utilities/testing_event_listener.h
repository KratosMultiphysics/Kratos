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
//

#pragma once

// System includes

// External includes
#include <gtest/gtest.h>

// Project includes

/*
 * ConfigurableEventListener provides a configurable event listener for the test output.
 * In Kratos this is used to remove the output from the tests
 * not executed by the main rank(0)
 * Inspiration from: From: https://gist.github.com/elliotchance/8215283
*/
class ConfigurableEventListener : public ::testing::TestEventListener
{
    protected:
        ::testing::TestEventListener * eventListener;

    public:

        /**
        * Swtiches to enable or disable the output of the test. Separated in different sections:
        * - showStart: Show the start of the test program.
        * - showIterations: Show the start of an iteration of the test program.
        * - showTestCases: Show the names of each test case.
        * - showTestNames: Show the names of each test.
        * - showSuccesses: Show each success.
        * - showInlineFailures: Show each failure as it occurs. You will also see it at the bottom after the full suite is run.
        * - showEnvironment: Show the setup of the global environment.
        * - showResult: Show the results of the test program.
        * - showEnd: Show the end of the test program.
        */

        /// Show the start of the test program.
        bool showStart;

        /// Show the start of an iteration of the test program.
        bool showIterations;

        /// Show the start of the test program.
        bool showTestCases;

        /// Show the start of the test program.
        bool showTestNames;

        /// Show the start of the test program.
        bool showSuccesses;

        /// Show the start of the test program.
        bool showInlineFailures;

        /// Show the start of the test program.
        bool showEnvironment;

        /// Show the start of the test program.
        bool showResult;

        /// Show the start of the test program.
        bool showEnd;

        explicit ConfigurableEventListener(::testing::TestEventListener* theEventListener) : eventListener(theEventListener)
        {
            showStart = true;
            showIterations = true;
            showTestCases = true;
            showTestNames = true;
            showSuccesses = true;
            showInlineFailures = true;
            showEnvironment = true;
            showResult = true;
            showEnd = true;
        }

        virtual ~ConfigurableEventListener()
        {
            delete eventListener;
        }

        virtual void OnTestProgramStart(const ::testing::UnitTest& unit_test) override
        {
            if(showStart) {
                eventListener->OnTestProgramStart(unit_test);
            }
        }

        virtual void OnTestIterationStart(const ::testing::UnitTest& unit_test, int iteration) override
        {
            if(showIterations) {
                eventListener->OnTestIterationStart(unit_test, iteration);
            }
        }

        virtual void OnEnvironmentsSetUpStart(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsSetUpStart(unit_test);
            }
        }

        virtual void OnEnvironmentsSetUpEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsSetUpEnd(unit_test);
            }
        }

        virtual void OnTestCaseStart(const ::testing::TestCase& test_case) override
        {
            if(showTestCases) {
                eventListener->OnTestCaseStart(test_case);
            }
        }

        virtual void OnTestStart(const ::testing::TestInfo& test_info) override
        {
            if(showTestNames) {
                eventListener->OnTestStart(test_info);
            }
        }

        virtual void OnTestPartResult(const ::testing::TestPartResult& result) override
        {
            if(showResult) {
                eventListener->OnTestPartResult(result);
            }
        }

        virtual void OnTestEnd(const ::testing::TestInfo& test_info) override
        {
            if((showInlineFailures && test_info.result()->Failed()) || (showSuccesses && !test_info.result()->Failed())) {
                eventListener->OnTestEnd(test_info);
            }
        }

        virtual void OnTestCaseEnd(const ::testing::TestCase& test_case) override
        {
            if(showTestCases) {
                eventListener->OnTestCaseEnd(test_case);
            }
        }

        virtual void OnEnvironmentsTearDownStart(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsTearDownStart(unit_test);
            }
        }

        virtual void OnEnvironmentsTearDownEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsTearDownEnd(unit_test);
            }
        }

        virtual void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration) override
        {
            if(showIterations) {
                eventListener->OnTestIterationEnd(unit_test, iteration);
            }
        }

        virtual void OnTestProgramEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnd) {
                eventListener->OnTestProgramEnd(unit_test);
            }
        }
};