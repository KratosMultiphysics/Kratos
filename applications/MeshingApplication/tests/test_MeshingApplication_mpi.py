# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
from test_mpi_metrics import TestMPIMetrics
from test_mpi_parmmg import TestMPIParMmg

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
    smallMPISuite.addTest(TestMPIMetrics('test_mpi_hessian'))
    if hasattr(KM.MeshingApplication,  "ParMmgProcess3D"):
        smallMPISuite.addTest(TestMPIParMmg('test_mpi_sphere'))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
