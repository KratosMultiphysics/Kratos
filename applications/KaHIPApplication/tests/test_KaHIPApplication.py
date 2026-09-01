# import Kratos
import KratosMultiphysics
import KratosMultiphysics.KaHIPApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the test classes to create the suites
from test_kahip_partitioner import TestKaHIPPartitioner
from test_kahip_modeler import TestKaHIPModeler

def AssembleTestSuites():
    """Populate the test suites to run.

    Returns
    -------
    suites : dict
        The set of suites with their test cases.
    """
    suites = KratosUnittest.KratosSuites

    # Small / fast tests
    small_suite = suites['small']
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        TestKaHIPPartitioner,
        TestKaHIPModeler,
    ]))

    # Nightly tests (include small)
    nightly_suite = suites['nightly']
    nightly_suite.addTests(small_suite)

    # All tests
    all_suite = suites['all']
    all_suite.addTests(nightly_suite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
