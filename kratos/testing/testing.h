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

class KRATOS_API(KRATOS_TEST_UTILS) KratosTestEnv : public ::testing::Environment 
{
    public:
        KratosTestEnv();
        ~KratosTestEnv() override {}
        void SetUp() override;
        void TearDown() override;
};

/*
 * Initializes the parallel testing environment. This is usefull for other tests depending on a parallel environment.
*/
class GTestMain {
    public:
        static int InitializeKernel(int argc, char* argv[]) {
            std::cout << "Kratos::Testing::GTestMain::InitializeKernel" << std::endl;
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
        // std::stringstream mStream;                                          // Stream to store the output of the tests and control visibility
        // std::streambuf * mCoutBuffer;
        // std::streambuf * mCerrBuffer;
        std::vector<KratosApplication::Pointer> mRegisteredApplications;    // List of applications loaded by the suit. TODO: Remove once every test includes its own suit
};

class KratosSensitivityTestSuite : public KratosCoreFastSuite {};
class KratosCoreGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreGeometryContainerFastSuite : public KratosCoreFastSuite {};
class KratosCoreNurbsGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosCoreCouplingGeometriesFastSuite : public KratosCoreFastSuite {};
class KratosExternalLibrariesFastSuite : public KratosCoreFastSuite {};
class KratosNonRectangularJacobianFastSuite : public KratosCoreFastSuite {};
class KratosCoreStressSuite : public KratosCoreFastSuite {};
class AdditiveManufacturingApplicationFastSuite: public KratosCoreFastSuite {};

// This classes are temporal and should be changed. Please see: GeoMechanicsApplication, StructuralMechanicsApplication or TrilinosApplication
// TODO: Remove once every test includes its own suit
class FluidDynamicsApplicationFastSuite : public KratosCoreFastSuite {};
class CompressiblePotentialApplicationFastSuite : public KratosCoreFastSuite {};
class KratosConstitutiveLawsFastSuite : public KratosCoreFastSuite {};
class KratosContactStructuralMechanicsFastSuite : public KratosCoreFastSuite {};
class KratosConvectionDiffusionFastSuite : public KratosCoreFastSuite {};
class KratosCSharpWrapperApplicationFastSuite : public KratosCoreFastSuite {};
class ExaquteSandboxApplicationFastSuite : public KratosCoreFastSuite {};
class FluidDynamicsBiomedicalApplicationFastSuite : public KratosCoreFastSuite {};
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

KRATOS_API(KRATOS_TEST_UTILS) DataCommunicator& GetDefaultDataCommunicator();

} // namespace Kratos::Testing
