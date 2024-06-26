//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kernel.h"
#include "tests/tests_utilities/test_suite.h"
#include "tests/tests_utilities/runtime_dependency_handler.h"

/*
 * This Fixture extends the KratosCoreFastSuite to include the runtime dependency handler.
 * This header is not ment to be used directly, but serve as a example on how to create a suite that handles runtime dependencies
 * in your own applications
 */
class KRATOS_API(KRATOS_TEST_UTILS) RuntimeDependencyExample: public KratosCoreFastSuite
{
    public:
        RuntimeDependencyExample(): KratosCoreFastSuite() {
            // The list of applications has to be a member, so that it has the same lifetime as the suite
            mRuntimeDependencyApplications = mRuntimeDependencies.CreateApplications();
            for (auto& r_dependency: mRuntimeDependencyApplications) {
                std::cout << "Registering " << r_dependency.first << std::endl;
                this->ImportApplicationIntoKernel(r_dependency.second);
            }
        }

        static void SetUpTestSuite() {
            // Dependencies should have a CreateApplication function to allow programatic import
            // The second argument is the library name, the extension (and the "lib" prefix, if needed) is added internally
            mRuntimeDependencies.LoadDependency("LinearSolversApplication", "KratosLinearSolversCore");
        }

        static void TearDownTestSuite() {
            mRuntimeDependencies.ReleaseDependencies();
        }

    private:

        std::unordered_map<std::string, KratosApplication::Pointer> mRuntimeDependencyApplications;

        static RuntimeDependencyHandler mRuntimeDependencies;
};