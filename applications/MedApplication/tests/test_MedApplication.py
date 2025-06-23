import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests of test-classes to create the suites
import test_med_model_part_io
import test_med_model_part_io_read_smp
import test_import_med_modeler


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

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites["small"]
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_med_model_part_io.TestMedModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_med_model_part_io_read_smp.TestMedModelPartIOReadSubModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_import_med_modeler.TestImportMedModeler]))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites["nightly"]
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites["all"]
    allSuite.addTests(nightSuite)

    return suites


if __name__ == "__main__":
    KratosUnittest.runTests(AssembleTestSuites())
