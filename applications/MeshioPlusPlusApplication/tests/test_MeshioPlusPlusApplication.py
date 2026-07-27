import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests of test-classes to create the suites
import test_meshioplusplus_io
import test_meshio_output_process
import test_meshio_input_modeler
import test_meshio_mesh_operations


def AssembleTestSuites():
    """Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    """

    suites = KratosUnittest.KratosSuites

    smallSuite = suites["small"]
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_meshioplusplus_io.TestMeshioPlusPlusIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_meshio_output_process.TestMeshioOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_meshio_input_modeler.TestMeshioInputModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_meshio_mesh_operations.TestMeshioPlusPlusMeshOperations]))

    nightSuite = suites["nightly"]
    nightSuite.addTests(smallSuite)

    allSuite = suites["all"]
    allSuite.addTests(nightSuite)

    return suites


if __name__ == "__main__":
    KratosUnittest.runTests(AssembleTestSuites())
