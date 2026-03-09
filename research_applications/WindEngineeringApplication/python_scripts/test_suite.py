# Kratos imports
import KratosMultiphysics.KratosUnittest as UnitTest

# STL imports
from enum import Enum


class SuiteFlags(Enum):
    ALL         = 0
    SMALL       = 1
    NIGHTLY     = 2
    VALIDATION  = 3
    MPI         = 4     # prepends suite names with "mpi_"
    MPI_ONLY    = 5     # same as above but won't get added to single process suites


class TestSuite(UnitTest.TestSuite):
    """Custom test suite class for sorting cases into suites automatically while globbing."""

    suite_flags = set([SuiteFlags.ALL])


    def addTest(self, test):
        if hasattr(test, "suite_flags") and (isinstance(test.suite_flags, set) or isinstance(test.suite_flags, list) or isinstance(test.suite_flags, tuple)):
            for flag in test.suite_flags:
                self.suite_flags.add(flag)
        UnitTest.TestSuite.addTest(self, test)