import os

import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as MetisApplication
    import KratosMultiphysics.TrilinosApplication as TrilinosApplication
except ImportError:
    raise Exception("KratosMPI could not be imported!")

# Import the tests or test_classes to create the suits
# process test_classes
from custom_process_tests import CustomProcessTest

# flow solver test_classes
from incompressible_potential_flow_solver_formulation_tests import IncompressiblePotentialFlowSolverFormulationTest
from monolithic_velocity_pressure_formulation_tests import MonolithicVelocityPressureFormulationTest
from fractional_step_velocity_pressure_formulation_tests import FractionalStepVelocityPressureFormulationTest

# turbulence model test_classes
### k-epsilon test_classes
from monolithic_k_epsilon_formulation_tests import MonolithicKEpsilonTest
from fractional_step_k_epsilon_formulation_tests import FractionalStepKEpsilonTest

### k_omega test_classes
from monolithic_k_omega_formulation_tests import MonolithicKOmegaTest
from fractional_step_k_omega_formulation_tests import FractionalStepKOmegaTest

### k_omega test_classes
from monolithic_k_omega_sst_formulation_tests import MonolithicKOmegaSSTTest
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

    # adding custom process tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CustomProcessTest]))

    # adding incompressible potential flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IncompressiblePotentialFlowSolverFormulationTest]))

    # adding monolithic flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicVelocityPressureFormulationTest]))

    # adding fractional step flow solver tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepVelocityPressureFormulationTest]))

    # adding monolithic k-epsilon high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKEpsilonTest]))
    # adding fractional step k-epsilon high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKEpsilonTest]))

    # adding monolithic k-omega high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaTest]))
    # adding fractional step k-omega high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaTest]))

    # adding monolithic k-omega-sst high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaSSTTest]))
    # adding fractional step k-omega-sst high re tests
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaSSTTest]))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites


if __name__ == '__main__':
    # this is required by the CI since, CI runs these tests from $KRATOS_HOME folder.
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
        KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
