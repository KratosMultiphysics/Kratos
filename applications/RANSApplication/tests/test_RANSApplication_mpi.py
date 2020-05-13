import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as MetisApplication
    import KratosMultiphysics.TrilinosApplication as TrilinosApplication
except ImportError:
    raise Exception("KratosMPI could not be imported!")

if KratosMultiphysics.ParallelEnvironment.GetDefaultSize() != 2:
    raise Exception("The MPI tests currently suport only being run with 2 processors!")

# Import the tests or test_classes to create the suits

# Shell tests
from evm_k_epsilon_tests import EvmKEpsilonTest
from incompressible_potential_flow_solver_formulation_tests import IncompressiblePotentialFlowSolverFormulationTest
from fractional_step_velocity_pressure_formulation_tests import FractionalStepVelocityPressureFormulationTest
from monolithic_velocity_pressure_formulation_tests import MonolithicVelocityPressureFormulationTest
from fractional_step_k_epsilon_high_re_formulation_tests import FractionalStepKEpsilonHighReTest
from monolithic_k_epsilon_high_re_formulation_tests import MonolithicKEpsilonHighReTest
from fractional_step_k_omega_formulation_tests import FractionalStepKOmegaTest
from monolithic_k_omega_formulation_tests import MonolithicKOmegaTest
from fractional_step_k_omega_sst_formulation_tests import FractionalStepKOmegaSSTTest
from monolithic_k_omega_sst_formulation_tests import MonolithicKOmegaSSTTest


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    ### Small MPI tests ########################################################
    # smallMPISuite = suites['mpi_small']

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTest((IncompressiblePotentialFlowSolverFormulationTest('testIncompressiblePotentialFlowMPI')))
    nightlyMPISuite.addTest((FractionalStepVelocityPressureFormulationTest('testFractionalStepVelocityPressureMPI')))
    nightlyMPISuite.addTest((MonolithicVelocityPressureFormulationTest('testMonolithicVelocityPressureMPI')))

    # ## fractional step tests
    nightlyMPISuite.addTest((FractionalStepKEpsilonHighReTest('testFractionalStepKEpsilonHighReAfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKEpsilonHighReTest('testFractionalStepKEpsilonHighReAfcVelocityMPI')))
    nightlyMPISuite.addTest((FractionalStepKEpsilonHighReTest('testFractionalStepKEpsilonHighReRfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKEpsilonHighReTest('testFractionalStepKEpsilonHighReRfcVelocityMPI')))

    nightlyMPISuite.addTest((FractionalStepKOmegaTest('testFractionalStepKOmegaAfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaTest('testFractionalStepKOmegaAfcVelocityMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaTest('testFractionalStepKOmegaRfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaTest('testFractionalStepKOmegaRfcVelocityMPI')))

    nightlyMPISuite.addTest((FractionalStepKOmegaSSTTest('testFractionalStepKOmegaSSTAfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaSSTTest('testFractionalStepKOmegaSSTAfcVelocityMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaSSTTest('testFractionalStepKOmegaSSTRfcTkeMPI')))
    nightlyMPISuite.addTest((FractionalStepKOmegaSSTTest('testFractionalStepKOmegaSSTRfcVelocityMPI')))

    ### monolithic tests
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcVelocityMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcVelocityMPI')))

    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcVelocityMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcVelocityMPI')))

    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTAfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTAfcVelocityMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcTkeMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcVelocityMPI')))


    nightlyMPISuite.addTest(EvmKEpsilonTest('testChannelFlowKEpsilonTransientMPI'))
    # nightlyMPISuite.addTest(EvmKEpsilonTest('testChannelFlowKEpsilonSteadyMPI'))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
