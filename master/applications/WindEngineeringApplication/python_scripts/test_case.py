# Kratos imports
import KratosMultiphysics.KratosUnittest as UnitTest

# Internal imports
from KratosMultiphysics.WindEngineeringApplication.test_suite import SuiteFlags, TestSuite


class TestCase(UnitTest.TestCase):
    """Custom test case class for sorting cases into suites automatically while globbing."""

    suite_flags = set([SuiteFlags.ALL])