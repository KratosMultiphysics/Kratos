# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import EmbeddedArtificialCompressibilityTest as TEmbeddedArtificialCompressibilityTest
from SmallTests import EmbeddedCouette2DTest as TEmbeddedCouette2DTest
from SmallTests import EmbeddedCouette3DTest as TEmbeddedCouette3DTest
from SmallTests import EmbeddedCouette2DImposedTest as TEmbeddedCouette2DImposedTest
from SmallTests import EmbeddedCouette3DImposedTest as TEmbeddedCouette3DImposedTest
from SmallTests import EmbeddedSlipReservoirTest as TEmbeddedSlipReservoirTest
from SmallTests import EmbeddedSlipBoundaryConditionTest as TEmbeddedSlipBoundaryConditionTest
from SmallTests import ManufacturedSolutionTest as TManufacturedSolutionTest
from SmallTests import NavierStokesWallConditionTest as TNavierStokesWallConditionTest

from buoyancy_test import BuoyancyTest
from volume_source_test import VolumeSourceTest
from fluid_element_test import FluidElementTest
from time_integrated_fluid_element_test import TimeIntegratedFluidElementTest
from darcy_channel_test import DarcyChannelTest
from embedded_reservoir_test import EmbeddedReservoirTest
from embedded_ausas_couette_test import EmbeddedAusasCouetteTest

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
    smallSuite.addTest(TEmbeddedCouette2DTest('test_execution'))
    smallSuite.addTest(TEmbeddedCouette2DImposedTest('test_execution'))
    smallSuite.addTest(TEmbeddedSlipBoundaryConditionTest('test_execution'))
    smallSuite.addTest(TManufacturedSolutionTest('test_execution'))
    smallSuite.addTest(TNavierStokesWallConditionTest('test_execution'))
    smallSuite.addTest(BuoyancyTest('testEulerian'))
    smallSuite.addTest(BuoyancyTest('testThermalExpansionCoefficient'))
    smallSuite.addTest(FluidElementTest('testCavity'))
    smallSuite.addTest(FluidElementTest('testCavityOSS'))
    smallSuite.addTest(TimeIntegratedFluidElementTest('testCavity'))
    smallSuite.addTest(TimeIntegratedFluidElementTest('testSymbolic'))
    smallSuite.addTest(DarcyChannelTest('testDarcyDensity'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir2D'))
    smallSuite.addTest(EmbeddedAusasCouetteTest('testEmbeddedAusasCouette2D'))
    smallSuite.addTest(EmbeddedAusasCouetteTest('testEmbeddedDevelopmentCouette2D'))
    #smallSuite.addTest(BuoyancyTest('testBFECC')) # I'm skipping this one, it varies too much between runs JC.

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(BuoyancyTest('validationEulerian'))
    #validationSuite.addTest(VolumeSourceTest('validationEulerian'))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    allSuite.addTest(DarcyChannelTest('testDarcyLinear'))
    allSuite.addTest(DarcyChannelTest('testDarcyNonLinear'))
    allSuite.addTest(TEmbeddedCouette3DTest('test_execution'))
    allSuite.addTest(TEmbeddedCouette3DImposedTest('test_execution'))
    allSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir3D'))
    allSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir3D'))
    allSuite.addTest(EmbeddedAusasCouetteTest('testEmbeddedAusasCouette3D'))
    smallSuite.addTest(EmbeddedAusasCouetteTest('testEmbeddedDevelopmentCouette3D'))

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
