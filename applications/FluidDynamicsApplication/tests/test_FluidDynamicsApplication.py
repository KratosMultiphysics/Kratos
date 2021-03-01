import subprocess
import os.path

# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites
from artificial_compressibility_test import ArtificialCompressibilityTest
from buoyancy_test import BuoyancyTest
from couette_flow_test import CouetteFlowTest
from darcy_channel_test import DarcyChannelTest
from embedded_piston_test import EmbeddedPistonTest
from embedded_reservoir_test import EmbeddedReservoirTest
from embedded_reservoir_discontinuous_test import EmbeddedReservoirDiscontinuousTest
from embedded_couette_flow_test import EmbeddedCouetteFlowTest
from embedded_couette_imposed_test import EmbeddedCouetteImposedTest
from embedded_velocity_inlet_emulation_test import EmbeddedVelocityInletEmulationTest
from fluid_element_test import FluidElementTest
from manufactured_solution_test import ManufacturedSolutionTest
from navier_stokes_wall_condition_test import NavierStokesWallConditionTest
from sod_shock_tube_test import SodShockTubeTest
from fluid_analysis_test import FluidAnalysisTest
from adjoint_fluid_test import AdjointFluidTest
from adjoint_vms_element_2d import AdjointVMSElement2D
from adjoint_vms_sensitivity_2d import AdjointVMSSensitivity2D
from hdf5_io_test import HDF5IOTest
from test_statistics_process import IntegrationPointStatisticsTest
from cfl_output_process_test import CFLOutputProcessTest
from test_flows_measuring_utility import FlowsMeasuringUtilityTest
from levelset_consistent_nodal_gradient_test import ConsistentLevelsetNodalGradientTest

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

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(CouetteFlowTest('testCouetteFlow2DSymbolicStokes'))
    smallSuite.addTest(CouetteFlowTest('testCouetteFlow2DWeaklyCompressibleNavierStokes'))
    smallSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokes2D'))
    smallSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokes2DSlip'))
    smallSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokesDiscontinuous2D'))
    smallSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed2D'))
    smallSuite.addTest(EmbeddedPistonTest('testEmbeddedPiston2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir2D'))
    smallSuite.addTest(EmbeddedVelocityInletEmulationTest('testEmbeddedVelocityInletEmulationEmbedded2D'))
    smallSuite.addTest(EmbeddedVelocityInletEmulationTest('testEmbeddedVelocityInletEmulationSymbolic2D'))
    smallSuite.addTest(NavierStokesWallConditionTest('testNavierStokesWallCondition'))
    smallSuite.addTest(FluidAnalysisTest('testSteadyAnalysisSmall'))
    #smallSuite.addTest(BuoyancyTest('testBFECC')) # I'm skipping this one, it varies too much between runs JC.

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(ArtificialCompressibilityTest('testArtificialCompressibility'))
    nightSuite.addTest(BuoyancyTest('testEulerian'))
    nightSuite.addTest(BuoyancyTest('testThermalExpansionCoefficient'))
    nightSuite.addTest(DarcyChannelTest('testDarcyDensity'))
    nightSuite.addTest(DarcyChannelTest('testDarcyLinear'))
    nightSuite.addTest(DarcyChannelTest('testDarcyNonLinear'))
    nightSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir3D'))
    nightSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir3D'))
    nightSuite.addTest(EmbeddedReservoirDiscontinuousTest('testEmbeddedReservoirDiscontinuous3D'))
    nightSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokes3D'))
    nightSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokes3DSlip'))
    nightSuite.addTest(EmbeddedCouetteFlowTest('testEmbeddedCouetteFlowEmbeddedWeaklyCompressibleNavierStokesDiscontinuous3D'))
    nightSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed3D'))
    nightSuite.addTest(FluidElementTest('testCavityQSASGS'))
    nightSuite.addTest(FluidElementTest('testCavityQSOSS'))
    nightSuite.addTest(FluidElementTest('testCavityDASGS'))
    nightSuite.addTest(FluidElementTest('testCavityDOSS'))
    nightSuite.addTest(FluidElementTest('testTimeIntegratedQSVMS'))
    nightSuite.addTest(FluidElementTest('testSymbolic'))
    nightSuite.addTest(FluidAnalysisTest('testFluidDynamicsAnalysis'))
    nightSuite.addTest(AdjointFluidTest('testCylinder'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSecondDerivativesLHS'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS1'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS2'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSensitivityMatrix'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testOneElement'))
    nightSuite.addTest(HDF5IOTest('testInputOutput'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCavity'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCylinder'))
    nightSuite.addTest(ConsistentLevelsetNodalGradientTest('testConsistentGradientSquare2D'))
    nightSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitASGS'))
    nightSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitASGSShockCapturing'))
    nightSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitOSS'))
    nightSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitOSSShockCapturing'))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IntegrationPointStatisticsTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CFLOutputProcessTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowsMeasuringUtilityTest]))


    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(BuoyancyTest('validationEulerian'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testCylinder'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testSteadyCylinder'))
    validationSuite.addTest(ManufacturedSolutionTest('testManufacturedSolution'))


    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
