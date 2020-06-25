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
from incompressible_potential_flow_solver_formulation_tests import IncompressiblePotentialFlowSolverFormulationTest
from monolithic_velocity_pressure_formulation_tests import MonolithicVelocityPressureFormulationTest
from monolithic_k_epsilon_high_re_formulation_tests import MonolithicKEpsilonHighReTest
from monolithic_k_omega_formulation_tests import MonolithicKOmegaTest
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
    nightlyMPISuite.addTest((MonolithicVelocityPressureFormulationTest('testMonolithicVelocityPressureSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicVelocityPressureFormulationTest('testMonolithicVelocityPressureTransientMPI')))

    # adding monolithic k-epsilon high re tests
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcTkeTransientMPI')))
    nightlyMPISuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcVelocityTransientMPI')))

    # adding monolithic k-omega tests
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcTkeTransientMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcVelocityTransientMPI')))

    # adding monolithic k-omega-sst tests
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTAfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTAfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcTkeSteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcVelocitySteadyMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcTkeTransientMPI')))
    nightlyMPISuite.addTest((MonolithicKOmegaSSTTest('testMonolithicKOmegaSSTRfcVelocityTransientMPI')))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
        KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
