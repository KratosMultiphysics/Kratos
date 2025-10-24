# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_perfect_plasticity_implementation_verification import TestPerfectPlasticityImplementationVerification

from test_factory import SimpleSmallDeformationPlasticityMCTest
from test_factory import SimpleSmallDeformationPlasticityVMTest
from test_factory import SimpleSmallDeformationPlasticityDPTest
from test_factory import SimpleSmallDeformationPlasticityTTest
from test_factory import BigCubeSmallDeformationPlasticityMCTest
from test_factory import BigCubeSmallDeformationPlasticityVMTest
from test_factory import BigCubeSmallDeformationPlasticityDPTest
from test_factory import BigCubeSmallDeformationPlasticityTTest
from test_factory import SerialParallelRuleOfMixturesCubeDamageTest
from test_factory import PlasticDamageTest
from test_factory import AnisotropyTest
from test_factory import Anisotropy2DTest
from test_factory import InitialStateInelasticityTest
from test_factory import InitialStateInelasticity2Test
from test_factory import SmallDeformationPlasticityTest
from test_factory import SimpleJ2PlasticityTest
from test_factory import TensileTestStructuralTest
from test_factory import HighCycleFatigueTest
from test_factory import PlasticDamageTest
from test_factory import AutomatedInitialDamageTest
from test_factory import TractionSeparationLawTest
from test_factory import CurveByPointsPlasticityTest
from test_factory import PlaneStressJ2Plasticity


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']

    smallSuite.addTest(SimpleSmallDeformationPlasticityVMTest('test_execution'))
    smallSuite.addTest(SimpleSmallDeformationPlasticityDPTest('test_execution'))
    smallSuite.addTest(SimpleSmallDeformationPlasticityTTest('test_execution'))
    smallSuite.addTest(BigCubeSmallDeformationPlasticityMCTest('test_execution'))
    smallSuite.addTest(BigCubeSmallDeformationPlasticityVMTest('test_execution'))
    smallSuite.addTest(BigCubeSmallDeformationPlasticityDPTest('test_execution'))
    smallSuite.addTest(BigCubeSmallDeformationPlasticityTTest('test_execution'))
    smallSuite.addTest(AnisotropyTest('test_execution'))
    smallSuite.addTest(Anisotropy2DTest('test_execution'))
    smallSuite.addTest(InitialStateInelasticityTest('test_execution'))
    smallSuite.addTest(InitialStateInelasticity2Test('test_execution'))
    smallSuite.addTest(SimpleJ2PlasticityTest('test_execution'))
    smallSuite.addTest(HighCycleFatigueTest('test_execution'))
    smallSuite.addTest(PlasticDamageTest('test_execution'))
    smallSuite.addTest(AutomatedInitialDamageTest('test_execution'))
    smallSuite.addTest(TractionSeparationLawTest('test_execution'))
    smallSuite.addTest(CurveByPointsPlasticityTest('test_execution'))
    smallSuite.addTest(PlaneStressJ2Plasticity('test_execution'))


    # Create a test suit with the selected tests (Nightly tests):
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPerfectPlasticityImplementationVerification]))
    nightSuite.addTest(SmallDeformationPlasticityTest('test_execution'))
    nightSuite.addTest(SimpleSmallDeformationPlasticityMCTest('test_execution'))
    nightSuite.addTest(SerialParallelRuleOfMixturesCubeDamageTest('test_execution'))
    nightSuite.addTest(PlasticDamageTest('test_execution'))

    ### Adding Validation Tests
    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(TensileTestStructuralTest('test_execution'))

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # Already contains the smallSuite
    validationSuite.addTests(allSuite) # Validation contains all

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
