import subprocess
import os.path


try:
    import sympy
    sympy_available = True
except:
    sympy_available = False
    print("Skipping tests that require sympy")

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
from embedded_couette_imposed_flow_test import EmbeddedCouetteImposedFlowTest
from embedded_velocity_inlet_emulation_test import EmbeddedVelocityInletEmulationTest
from fluid_element_test import FluidElementTest
from manufactured_solution_test import ManufacturedSolutionTest
from navier_stokes_wall_condition_test import NavierStokesWallConditionTest
from sod_shock_tube_test import SodShockTubeTest
from fluid_analysis_test import FluidAnalysisTest
from adjoint_fluid_test import AdjointFluidTest
from adjoint_vms_element_2d import AdjointVMSElement2D
from adjoint_vms_sensitivity_2d import AdjointVMSSensitivity2D
from adjoint_qsvms_sensitivity_2d import AdjointQSVMSSensitivity2D
from hdf5_io_test import HDF5IOTest
from test_statistics_process import IntegrationPointStatisticsTest
from test_flows_measuring_utility import FlowsMeasuringUtilityTest
from levelset_consistent_nodal_gradient_test import ConsistentLevelsetNodalGradientTest
from adjoint_conditions import TestAdjointMonolithicWallCondition
from test_fluid_auxiliary_utilities import FluidAuxiliaryUtilitiesTest
from test_navier_stokes_compressible_explicit_solver import NavierStokesCompressibleExplicitSolverTest
from two_fluid_mass_conservation_source_test import TwoFluidMassConservationTest
from apply_compressible_navier_stokes_boundary_conditions_process_test import ApplyMachDependentBoundaryConditionsTest
if sympy_available:
    from compressible_navier_stokes_symbolic_generator_formulation_test import CompressibleNavierStokesSymbolicGeneratorFormulationTest
from compressible_slip_wall_process_test import TestCompressibleSlipWallProcess
from compute_pressure_coefficient_process_test import ComputePressureCoefficientProcessTest
from compute_drag_process_test import ComputeDragProcessTest
from test_compute_y_plus_process import ComputeYPlusProcessTest
from test_fluid_computation_processes import FluidComputationProcessesTest
from slip_spurious_tangential_correction_test import SlipSpuriousTangentialCorrectionTest
from apply_wall_law_process_test import ApplyWallLawProcessTest

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
    smallSuite.addTest(EmbeddedCouetteImposedFlowTest('testEmbeddedCouetteImposedFlowEmbeddedWeaklyCompressibleNavierStokes2D'))
    smallSuite.addTest(EmbeddedPistonTest('testEmbeddedPiston2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir2D'))
    smallSuite.addTest(NavierStokesWallConditionTest('testNavierStokesWallCondition'))
    smallSuite.addTest(FluidAnalysisTest('testSteadyAnalysisSmall'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([EmbeddedVelocityInletEmulationTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointMonolithicWallCondition]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointMonolithicWallCondition]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ApplyMachDependentBoundaryConditionsTest]))
    #smallSuite.addTest(BuoyancyTest('testBFECC')) # I'm skipping this one, it varies too much between runs JC.
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCompressibleSlipWallProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ComputePressureCoefficientProcessTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ComputeDragProcessTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ComputeYPlusProcessTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SlipSpuriousTangentialCorrectionTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ApplyWallLawProcessTest]))

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
    nightSuite.addTest(EmbeddedCouetteImposedFlowTest('testEmbeddedCouetteImposedFlowEmbeddedWeaklyCompressibleNavierStokes3D'))
    nightSuite.addTest(FluidAnalysisTest('testFluidDynamicsAnalysis'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSecondDerivativesLHS'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS1'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS2'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSensitivityMatrix'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testOneElement'))
    nightSuite.addTest(AdjointQSVMSSensitivity2D('testOneElement'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testTwoElementsSlipSteady'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testTwoElementsSlipBossak'))
    nightSuite.addTest(HDF5IOTest('testInputOutput'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCavity'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCylinder'))
    nightSuite.addTest(ConsistentLevelsetNodalGradientTest('testConsistentGradientSquare2D'))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FluidElementTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IntegrationPointStatisticsTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowsMeasuringUtilityTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FluidAuxiliaryUtilitiesTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TwoFluidMassConservationTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([NavierStokesCompressibleExplicitSolverTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FluidComputationProcessesTest]))

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(BuoyancyTest('validationEulerian'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testCylinder'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testSteadyCylinder'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testSlipNormCylinder'))
    validationSuite.addTest(AdjointVMSSensitivity2D('testSlipSteadyNormCylinder'))
    validationSuite.addTest(AdjointQSVMSSensitivity2D('testCylinder'))
    validationSuite.addTest(AdjointQSVMSSensitivity2D('testSteadyCylinder'))
    validationSuite.addTest(ManufacturedSolutionTest('testManufacturedSolution'))
    #FIXME: MOVE BACK THE SOD TO NIGHT ONCE WE FIX THE NIGHTLY BUILD ISSUE
    validationSuite.addTest(AdjointFluidTest('testCylinder'))
    validationSuite.addTest(AdjointFluidTest('testSlipCylinder'))
    validationSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitASGS'))
    validationSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitASGSShockCapturing'))
    validationSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitOSS'))
    validationSuite.addTest(SodShockTubeTest('testSodShockTubeExplicitOSSShockCapturing'))
    if sympy_available:
        validationSuite.addTest(CompressibleNavierStokesSymbolicGeneratorFormulationTest('testSymbolicTriangle'))
        validationSuite.addTest(CompressibleNavierStokesSymbolicGeneratorFormulationTest('testSymbolicTetrahedron'))
        validationSuite.addTest(CompressibleNavierStokesSymbolicGeneratorFormulationTest('testSymbolicQuadrilateral'))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
        KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
