# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import EmbeddedArtificialCompressibilityTest as TEmbeddedArtificialCompressibilityTest

from buoyancy_test import BuoyancyTest
from darcy_channel_test import DarcyChannelTest
from embedded_reservoir_test import EmbeddedReservoirTest
from embedded_couette_test import EmbeddedCouetteTest
from embedded_couette_imposed_test import EmbeddedCouetteImposedTest
from fluid_element_test import FluidElementTest
from manufactured_solution_test import ManufacturedSolutionTest
from navier_stokes_wall_condition_test import NavierStokesWallConditionTest
from time_integrated_fluid_element_test import TimeIntegratedFluidElementTest
from volume_source_test import VolumeSourceTest

## NIGTHLY TESTS

## VALIDATION TESTS

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(TEmbeddedArtificialCompressibilityTest('test_execution'))
    smallSuite.addTest(NavierStokesWallConditionTest('testNavierStokesWallCondition'))
    smallSuite.addTest(ManufacturedSolutionTest('testManufacturedSolution'))
    smallSuite.addTest(BuoyancyTest('testEulerian'))
    smallSuite.addTest(BuoyancyTest('testThermalExpansionCoefficient'))
    smallSuite.addTest(FluidElementTest('testCavity'))
    smallSuite.addTest(FluidElementTest('testCavityOSS'))
    smallSuite.addTest(TimeIntegratedFluidElementTest('testCavity'))
    smallSuite.addTest(TimeIntegratedFluidElementTest('testSymbolic'))
    smallSuite.addTest(DarcyChannelTest('testDarcyDensity'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedDevelopmentCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipDevelopmentCouette2D'))
    smallSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed2D'))
    #smallSuite.addTest(BuoyancyTest('testBFECC')) # I'm skipping this one, it varies too much between runs JC.

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(BuoyancyTest('validationEulerian'))
    validationSuite.addTest(VolumeSourceTest('validationEulerian'))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    allSuite.addTest(DarcyChannelTest('testDarcyLinear'))
    allSuite.addTest(DarcyChannelTest('testDarcyNonLinear'))
    allSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir3D'))
    allSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir3D'))
    allSuite.addTest(EmbeddedCouetteTest('testEmbeddedCouette3D'))
    allSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipCouette3D'))
    allSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasCouette3D'))
    allSuite.addTest(EmbeddedCouetteTest('testEmbeddedDevelopmentCouette3D'))
    allSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipDevelopmentCouette3D'))
    allSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed3D'))

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
