# Kratos imports
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.WindEngineeringApplication.test_suite import SuiteFlags, TestSuite

# STL imports
import pathlib


class TestLoader(UnitTest.TestLoader):
    @property
    def suiteClass(self):
        return TestSuite


def AssembleTestSuites(enable_mpi=False):
    """ Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    """

    static_suites = UnitTest.KratosSuites

    # Test cases will be organized into lists first, then loaded into their
    # corresponding suites all at once
    local_cases = {}
    for key in static_suites.keys():
        local_cases[key] = []

    # Glob all test cases in this application
    this_directory = pathlib.Path(__file__).absolute().parent
    test_loader = TestLoader()
    all_tests = test_loader.discover(this_directory)

    # Sort globbed test cases into lists based on their suite flags
    #   flags correspond to entries in KratosUnittest.TestSuites
    #   (small, nightly, all, validation)
    #
    #   Cases with the 'mpi' flag are added to mpi suites as well as their corresponding normal suites.
    #   Cases with the 'mpi_only' flag are not added to normal suites.
    for test_case in all_tests:
        suite_flags = set(test_case.suite_flags)

        # Check whether the test case has a flag for mpi
        mpi = SuiteFlags.MPI in suite_flags
        mpi_only = SuiteFlags.MPI_ONLY in suite_flags

        # Don't add the test if its mpi-exclusive and mpi is not enabled
        if (not enable_mpi) and mpi_only:
            continue

        # Remove mpi flags
        if mpi:
            suite_flags.remove(SuiteFlags.MPI)

        if mpi_only:
            suite_flags.remove(SuiteFlags.MPI_ONLY)

        # Add case to the corresponding suites
        for suite_flag in suite_flags:
            local_cases[suite_flag.name.lower()].append(test_case)
            if mpi or mpi_only:
                local_cases["mpi_" + suite_flag.name.lower()].append(test_case)

        # Put test in 'all' if it isn't already there
        if not (SuiteFlags.ALL in suite_flags):
            if not mpi_only:
                local_cases["all"].append(test_case)
            if mpi or mpi_only:
                local_cases["mpi_all"].append(test_case)

    # Load all sorted cases into the global suites
    for suite_name, test_cases in local_cases.items():
        static_suites[suite_name].addTests(test_cases)

    return static_suites


def Run(enable_mpi=False):
    UnitTest.runTests(AssembleTestSuites(enable_mpi=enable_mpi))


if __name__ == "__main__":
    Run(enable_mpi=False)