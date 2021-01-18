# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

if KratosMultiphysics.ParallelEnvironment.GetDefaultSize() != 2:
    raise Exception("The MPI tests currently support only being run with 2 processors!")

# Import the tests or test_classes to create the suits

# Shell tests
from structural_mechanics_test_factory_mpi import ShellT3AndQ4LinearStaticStructPinchedCylinderTests as TShellT3AndQ4LinearStaticStructPinchedCylinderTests


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
    smallMPISuite = suites['mpi_small']

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTest(TShellT3AndQ4LinearStaticStructPinchedCylinderTests('test_execution'))
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
