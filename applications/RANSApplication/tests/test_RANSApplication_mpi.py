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

from fractional_step_velocity_pressure_formulation_tests import FractionalStepVelocityPressureFormulationTest
from fractional_step_k_epsilon_high_re_formulation_tests import FractionalStepKEpsilonHighReTest
from fractional_step_k_omega_formulation_tests import FractionalStepKOmegaTest
from fractional_step_k_omega_sst_formulation_tests import FractionalStepKOmegaSSTTest

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

        # addign incompressible potential flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IncompressiblePotentialFlowSolverFormulationTest]))

    # adding monolithic flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicVelocityPressureFormulationTest]))

    # adding monolithic k-epsilon high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKEpsilonHighReTest]))

    # adding monolithic k-omega tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaTest]))

    # adding monolithic k-omega-sst tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaSSTTest]))

    # adding fractional step flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepVelocityPressureFormulationTest]))

    # adding fractional step k-epsilon high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKEpsilonHighReTest]))

    # adding fractional step k-omega tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaTest]))

    # adding fractional step k-omega-sst tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaSSTTest]))

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
