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
//                   Jordi Cotela Dalmau
//
//

#pragma once

// System includes

// External includes
#include <gtest/gtest.h>

// Project includes
#include "includes/kernel.h"

/*
 * This Fixture creates a new kernel instance for kratos, so the test is able to interact with the database.
 * Its called this way to that all tests belong to a existing kernel fixture
*/
class KRATOS_API(KRATOS_TEST_UTILS) KratosCoreFastSuite : public ::testing::Test
{
    public:
        void SetUp() override;
        void TearDown() override;

    protected:
        KratosCoreFastSuite(): mKernel() {
            for (auto && appInitializer: mApplicationInitializerList) {
                appInitializer(mRegisteredApplications, mKernel);
            }
        }

        ~KratosCoreFastSuite() {}

        void ImportApplicationIntoKernel(KratosApplication::Pointer pNewApplication){
            if (!mKernel.IsImported(pNewApplication->Name())) {
                mKernel.ImportApplication(pNewApplication);
            }
        }

    private:
        Kratos::Kernel mKernel;
        // std::stringstream mStream;                                       // Stream to store the output of the tests and control visibility
        // std::streambuf * mCoutBuffer;
        // std::streambuf * mCerrBuffer;
        std::vector<KratosApplication::Pointer> mRegisteredApplications;    // List of applications loaded by the suit. TODO: Remove once every test includes its own suit
};