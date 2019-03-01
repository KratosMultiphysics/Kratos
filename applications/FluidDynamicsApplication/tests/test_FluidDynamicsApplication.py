# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

import subprocess

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites
from artificial_compressibility_test import ArtificialCompressibilityTest
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
from fluid_analysis_test import FluidAnalysisTest
from adjoint_fluid_test import AdjointFluidTest
from adjoint_vms_element_2d import AdjointVMSElement2D
from adjoint_vms_sensitivity_2d import AdjointVMSSensitivity2D
from hdf5_io_test import HDF5IOTest
from test_statistics_process import IntegrationPointStatisticsTest
from test_ale_fluid_solver import ALEFluidSolverTest

import run_cpp_unit_tests

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
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedDevelopmentCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipDevelopmentCouette2D'))
    smallSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasDevelopmentCouette2D'))
    smallSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedReservoir2D'))
    smallSuite.addTest(EmbeddedReservoirTest('testEmbeddedSlipReservoir2D'))
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
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedCouette3D'))
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipCouette3D'))
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasCouette3D'))
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedDevelopmentCouette3D'))
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedSlipDevelopmentCouette3D'))
    nightSuite.addTest(EmbeddedCouetteTest('testEmbeddedAusasDevelopmentCouette3D'))
    nightSuite.addTest(EmbeddedCouetteImposedTest('testEmbeddedCouetteImposed3D'))
    nightSuite.addTest(FluidElementTest('testCavityQSASGS'))
    nightSuite.addTest(FluidElementTest('testCavityQSOSS'))
    nightSuite.addTest(FluidElementTest('testCavityDASGS'))
    nightSuite.addTest(FluidElementTest('testCavityDOSS'))
    nightSuite.addTest(ManufacturedSolutionTest('testManufacturedSolution'))
    nightSuite.addTest(TimeIntegratedFluidElementTest('testCavity'))
    nightSuite.addTest(TimeIntegratedFluidElementTest('testSymbolic'))
    nightSuite.addTest(FluidAnalysisTest('testFluidDynamicsAnalysis'))
    nightSuite.addTest(AdjointFluidTest('testCylinder'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSecondDerivativesLHS'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS1'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateFirstDerivativesLHS2'))
    nightSuite.addTest(AdjointVMSElement2D('testCalculateSensitivityMatrix'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testCylinder'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testOneElement'))
    nightSuite.addTest(AdjointVMSSensitivity2D('testSteadyCylinder'))
    nightSuite.addTest(HDF5IOTest('testInputOutput'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCavity'))
    nightSuite.addTest(FluidAnalysisTest('testSteadyCylinder'))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IntegrationPointStatisticsTest]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ALEFluidSolverTest]))


    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(BuoyancyTest('validationEulerian'))
    validationSuite.addTest(VolumeSourceTest('validationEulerian'))


    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
    try:
        import KratosMultiphysics.mpi as KratosMPI
        import KratosMultiphysics.MetisApplication as MetisApplication
        import KratosMultiphysics.TrilinosApplication as TrilinosApplication
        p = subprocess.Popen(["mpiexec", "-np", "2", "python3", "test_FluidDynamicsApplication_mpi.py"], stdout=subprocess.PIPE)
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    except ImportError:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "mpi is not available!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")