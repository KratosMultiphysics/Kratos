# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import from Test Factories (with general analysis flows)
from particle_mechanics_test_factory import BeamCantileverStaticLinearElasticPointLoad2DTriTest as TBeamCantileverStaticLinearElasticPointLoad2DTriTest
from particle_mechanics_test_factory import BeamCantileverStaticLinearElasticLineLoad2DQuadTest as TBeamCantileverStaticLinearElasticLineLoad2DQuadTest
from particle_mechanics_test_factory import BeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest as TBeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest
from particle_mechanics_test_factory import BeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest as TBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest
from particle_mechanics_test_factory import CooksMembraneCompressibleTest as TCooksMembraneCompressibleTest
from particle_mechanics_test_factory import CooksMembraneUPCompressibleTest as TCooksMembraneUPCompressibleTest
from particle_mechanics_test_factory import CooksMembraneUPIncompressibleTest as TCooksMembraneUPIncompressibleTest
from particle_mechanics_test_factory import CLLinearElastic3DQuadTest as TCLLinearElastic3DQuadTest

# Import from Test Factories (with different analysis flows)
from test_static_loading_conditions_point import TestStaticLoadingConditionsPoint as TTestStaticLoadingConditionsPoint
from test_static_loading_conditions_line  import TestStaticLoadingConditionsLine  as TTestStaticLoadingConditionsLine


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

    ### Create a test suit with the selected tests (Small tests):
    ### These tests have to be very fast!
    ### Execution time << 1 sec on a regular PC !!!

    ## These tests are executed by the continuous integration tool
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestStaticLoadingConditionsPoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestStaticLoadingConditionsLine]))

    ## These tests are executed in the nightly build
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    nightSuite.addTest(TBeamCantileverStaticLinearElasticPointLoad2DTriTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticLineLoad2DQuadTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest('test_execution'))

    nightSuite.addTest(TCooksMembraneCompressibleTest('test_execution'))
    nightSuite.addTest(TCooksMembraneUPCompressibleTest('test_execution'))
    nightSuite.addTest(TCooksMembraneUPIncompressibleTest('test_execution'))
    nightSuite.addTest(TCLLinearElastic3DQuadTest('test_execution'))

    ### Adding Validation Tests
    ## For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    ### Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
