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

namespace Kratos::Testing 
{

KRATOS_API(KRATOS_TEST_UTILS) extern std::vector<std::function<void(std::vector<KratosApplication::Pointer> &, Kratos::Kernel &)>> mApplicationInitializerList;

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

// Define some suits that are needed by core tests
class KratosSensitivityTestSuite : public KratosCoreFastSuite {};
class KratosCoreGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreGeometryContainerFastSuite : public KratosCoreFastSuite {};
class KratosCoreNurbsGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreCouplingGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosExternalLibrariesFastSuite : public KratosCoreFastSuite {};
class KratosNonRectangularJacobianFastSuite : public KratosCoreFastSuite {};
class KratosCoreStressSuite : public KratosCoreFastSuite {};

////////////////////////////////////////////////////////////////////////////////////////
//// ALL THIS CLASSES NEED TO BE REMOVED OR MOVED TO CORE_TEST_SUTES.H ONCE MERGED  ////
////////////////////////////////////////////////////////////////////////////////////////

// This classes are temporal and should be changed. Please see: GeoMechanicsApplication, StructuralMechanicsApplication or TrilinosApplication
// TODO: Remove once every test includes its own suit
class CompressiblePotentialApplicationFastSuite : public KratosCoreFastSuite {};
class KratosConstitutiveLawsFastSuite : public KratosCoreFastSuite {};
class KratosContactStructuralMechanicsFastSuite : public KratosCoreFastSuite {};
class KratosConvectionDiffusionFastSuite : public KratosCoreFastSuite {};
class KratosCSharpWrapperApplicationFastSuite : public KratosCoreFastSuite {};
class ExaquteSandboxApplicationFastSuite : public KratosCoreFastSuite {};
class FSIApplicationFastSuite : public KratosCoreFastSuite {};
class KratosIgaFastSuite : public KratosCoreFastSuite {};
class KratosIgaFast5PSuite : public KratosCoreFastSuite {};
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
class ShallowWaterApplicationFastSuite : public KratosCoreFastSuite {};

// TODO: those should be in its own mpi suite
class KratosMappingApplicationMPITestSuite : public KratosCoreFastSuite {};

} // namespace Kratos::Testing