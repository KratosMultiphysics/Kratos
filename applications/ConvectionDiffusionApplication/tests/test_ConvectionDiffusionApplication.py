# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

# Import the tests or test_classes to create the suites

##### SELF-CONTAINED TESTS #####
from bfecc_convection_test import BFECCConvectionTest
from source_term_test import SourceTermTest
from thermal_coupling_test import ThermalCouplingTest
from test_apply_thermal_face_process import ApplyThermalFaceProcessTest
from adjoint_heat_diffusion_test import AdjointHeatDiffusionTest

##### SMALL TESTS #####
from convection_diffusion_test_factory import BasicConvectionDiffusionStationaryTest as TBasicConvectionDiffusionStationaryTest
from convection_diffusion_test_factory import BasicConvectionDiffusionTransientTest as TBasicConvectionDiffusionTransientTest
from convection_diffusion_test_factory import BasicDiffusionStationaryTest as TBasicDiffusionStationaryTest
from convection_diffusion_test_factory import SimpleThermoMechanicalTest as TSimpleThermoMechanicalTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # These tests are executed by the continuous integration tool, so they have to be very fast!
    # Execution time << 1 sec on a regular PC !!!
    # If the tests in the smallSuite take too long then merging to master will not be possible!
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    nightSuite = suites['nightly'] # These tests are executed in the nightly build

    ### Adding the self-contained tests
    smallSuite.addTest(BFECCConvectionTest('testBFECCConvection'))
    smallSuite.addTest(BFECCConvectionTest('testBFECCElementalLimiterConvection'))
    smallSuite.addTest(SourceTermTest('testPureDiffusion'))
    smallSuite.addTest(SourceTermTest('testDiffusionDominated'))
    smallSuite.addTest(SourceTermTest('testConvectionDominated'))
    smallSuite.addTest(SourceTermTest('testReaction'))
    smallSuite.addTest(ThermalCouplingTest('testDirichletNeumann'))
    smallSuite.addTest(ApplyThermalFaceProcessTest('testThermalFaceProcess'))
    smallSuite.addTest(AdjointHeatDiffusionTest('testAdjointHeatDiffusion'))
    smallSuite.addTest(AdjointHeatDiffusionTest('testAdjointHeatDiffusionWithSourceTerm'))

    ### Adding Small Tests
    smallSuite.addTest(TBasicConvectionDiffusionStationaryTest('test_execution'))
    smallSuite.addTest(TBasicConvectionDiffusionTransientTest('test_execution'))
    smallSuite.addTest(TBasicDiffusionStationaryTest('test_execution'))
    smallSuite.addTest(TSimpleThermoMechanicalTest('test_execution'))

    # Create a test suite with the selected tests plus all small tests
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
