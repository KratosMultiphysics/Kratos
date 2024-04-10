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
//
//

#pragma once

// System includes

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "includes/kernel.h"
#include "includes/expect.h"                    // Includes the expects from gtest and gmock adapted to kratos checks.
#include "includes/data_communicator.h"
#include "testing/test_skipped_exception.h"     // Macros and exception class used to skip tests.

#define KRATOS_TEST_CASE(A) TEST_F(KratosCoreFastSuite, A)
#define KRATOS_TEST_CASE_IN_SUITE(A, B) TEST_F(B, A)
#define KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(A, B) TEST_F(B, A)

namespace Kratos::Testing {

/*
 * This Fixture creates a new kernel instance for kratos, so the test is able to interact with the database.
 * Its called this way to that all tests belong to a existing kernel fixture
*/
class KratosCoreFastSuite : public ::testing::Test 
{
    protected:
        KratosCoreFastSuite(): mKernel() {}
        ~KratosCoreFastSuite() {}

        Kratos::Kernel mKernel;
};

class KratosSensitivityTestSuite : public KratosCoreFastSuite {};
class KratosStructuralMechanicsFastSuite: public KratosCoreFastSuite {};
class KratosCoreGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreGeometryContainerFastSuite : public KratosCoreFastSuite {};
class KratosCoreNurbsGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreCouplingGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosExternalLibrariesFastSuite : public KratosCoreFastSuite {};
class KratosNonRectangularJacobianFastSuite : public KratosCoreFastSuite {};
class MeshMovingApplicationFastSuite : public KratosCoreFastSuite {};
class KratosStatisticsFastSuite: public KratosCoreFastSuite {};
class KratosCoreStressSuite : public KratosCoreFastSuite {};

// TODO: To be removed
class FluidDynamicsApplicationFastSuite : public KratosCoreFastSuite {};

KRATOS_API(KRATOS_TEST_UTILS) DataCommunicator& GetDefaultDataCommunicator();

}
