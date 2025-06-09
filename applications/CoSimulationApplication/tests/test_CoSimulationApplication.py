# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

from co_simulation_test_factory import TestTinyFetiCoSimulationCases
from co_simulation_test_factory import TestSmallCoSimulationCases
from co_simulation_test_factory import TestCoSimulationCases
from test_function_callback_utility import TestGenericCallFunction
from test_ping_pong_coupling import TestPingPong
from test_processes import TestCreatePointBasedEntitiesProcess
from test_mok_fsi import TestMokFSI
from test_coupling_interface_data import TestCouplingInterfaceData
from test_data_transfer_operators import TestDataTransferOperators
from test_coupling_operations import TestScalingOperation
from test_flower_coupling import TestFLOWerCoupling
from test_sdof_solver import TestSdofSolver
from test_sdof_static_solver import TestSdofStaticSolver
from test_convergence_criteria import TestConvergenceCriteria
from test_convergence_criteria import TestConvergenceCriteriaWrapper
from test_convergence_accelerators import TestConvergenceAcceleratorWrapper
from test_co_simulation_coupled_solver import TestCoupledSolverGetSolver
from test_co_simulation_coupled_solver import TestCoupledSolverModelAccess
from test_co_simulation_coupled_solver import TestCoupledSolverPassingModel
from test_co_simulation_coupled_solver import TestCoupledSolverCouplingInterfaceDataAccess
from test_model_part_utilties import TestModelPartUtiliites
from test_thermal_rom_co_sim import TestThermalRomCoSim
from test_3d_1d_data_transfer_process import Test3D1DDataTransferProcessBlock
from test_3d_1d_data_transfer_process import Test3D1DDataTransferProcessTorus

from test_cosim_EMPIRE_API import TestCoSim_EMPIRE_API
from test_co_sim_io_py_exposure import TestCoSimIOPyExposure
from test_co_sim_io_py_exposure import TestCoSimIOPyExposure_aux_tests
from test_kratos_co_sim_io import TestKratosCoSimIO


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
    ################################################################################
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGenericCallFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCreatePointBasedEntitiesProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSim_EMPIRE_API]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimIOPyExposure_aux_tests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCouplingInterfaceData]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDataTransferOperators]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestScalingOperation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSdofSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSdofStaticSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConvergenceCriteria]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConvergenceCriteriaWrapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoupledSolverGetSolver]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoupledSolverModelAccess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoupledSolverPassingModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoupledSolverCouplingInterfaceDataAccess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModelPartUtiliites]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPingPong]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConvergenceAcceleratorWrapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTinyFetiCoSimulationCases]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestThermalRomCoSim]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([Test3D1DDataTransferProcessBlock]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([Test3D1DDataTransferProcessTorus]))


    ################################################################################
    nightSuite = suites['nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSmallCoSimulationCases]))
    
    # This one has errors in GCC
    nightSuite.addTest(TestMokFSI('test_mok_fsi_mvqn'))
    nightSuite.addTest(TestMokFSI('test_mok_fsi_aitken'))
    
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimIOPyExposure]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestKratosCoSimIO]))

    nightSuite.addTests(smallSuite)

    ################################################################################
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimulationCases]))
    validationSuite.addTest(TestMokFSI('test_mok_fsi_mvqn_external_structure'))
    validationSuite.addTest(TestMokFSI('test_mok_fsi_mvqn_external_structure_remote_controlled'))
    # validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFLOWerCoupling]))

    ################################################################################
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
