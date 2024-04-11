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

KRATOS_API(KRATOS_TEST_UTILS) extern std::vector<std::function<void(std::vector<KratosApplication::Pointer> &, Kratos::Kernel &)>> mApplicationInitializerList;

/*
 * This Fixture creates a new kernel instance for kratos, so the test is able to interact with the database.
 * Its called this way to that all tests belong to a existing kernel fixture
*/
class KratosCoreFastSuite : public ::testing::Test 
{
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
        std::vector<KratosApplication::Pointer> mRegisteredApplications; // This needs to be declared <!<!AFTER!>!> the kernel
};

class KratosSensitivityTestSuite : public KratosCoreFastSuite {};
class KratosCoreGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreGeometryContainerFastSuite : public KratosCoreFastSuite {};
class KratosCoreNurbsGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreCouplingGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosExternalLibrariesFastSuite : public KratosCoreFastSuite {};
class KratosNonRectangularJacobianFastSuite : public KratosCoreFastSuite {};
class KratosCoreStressSuite : public KratosCoreFastSuite {};

// This classes are temporal and should be changed. Please see: GeoMechanicsApplication, StructuralMechanicsApplication or TrilinosApplication
class FluidDynamicsApplicationFastSuite : public KratosCoreFastSuite {};
class CompressiblePotentialApplicationFastSuite : public KratosCoreFastSuite {};
class KratosConstitutiveLawsFastSuite : public KratosCoreFastSuite {};
class KratosContactStructuralMechanicsFastSuite : public KratosCoreFastSuite {};
class KratosConvectionDiffusionFastSuite : public KratosCoreFastSuite {};
class KratosCosimulationFastSuite : public KratosCoreFastSuite {};
class KratosCSharpWrapperApplicationFastSuite : public KratosCoreFastSuite {};
class ExaquteSandboxApplicationFastSuite : public KratosCoreFastSuite {};
class FluidDynamicsBiomedicalApplicationFastSuite : public KratosCoreFastSuite {};
class FSIApplicationFastSuite : public KratosCoreFastSuite {};
class KratosIgaFastSuite : public KratosCoreFastSuite {};
class KratosMappingApplicationSerialTestSuite : public KratosCoreFastSuite {};
class KratosMedFastSuite : public KratosCoreFastSuite {};
class MeshMovingApplicationFastSuite : public KratosCoreFastSuite {};
class KratosMeshingApplicationFastSuite : public KratosCoreFastSuite {};
class KratosRansFastSuite : public KratosCoreFastSuite {};
class RomApplicationFastSuite : public KratosCoreFastSuite {};
class KratosSolidMechanicsFastSuite : public KratosCoreFastSuite {};
class KratosStatisticsFastSuite : public KratosCoreFastSuite {};
class KratosWindEngineeringFastSuite : public KratosCoreFastSuite {};
class KratosMPMFastSuite : public KratosCoreFastSuite {};
class KratosHDF5TestSuite : public KratosCoreFastSuite {};
class KratosMappingApplicationMPITestSuite : public KratosCoreFastSuite {};

KRATOS_API(KRATOS_TEST_UTILS) DataCommunicator& GetDefaultDataCommunicator();

}
