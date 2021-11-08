# import Kratos

import os
import sys

sys.path.append(os.path.join('..', '..', '..'))

from KratosMultiphysics import *
from KratosMultiphysics.GeoMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosGeoMechanicsGeneralTests
from test_excavation import KratosGeoMechanicsExcavationTests
from test_interface import KratosGeoMechanicsInterfaceTests
from test_reset_displacement import KratosGeoMechanicsResetDisplacementTests
from test_benchmark_set_1 import KratosGeoMechanicsBenchmarkSet1
from test_benchmark_set_2 import KratosGeoMechanicsBenchmarkSet2
from test_soil_structure_interactions import KratosGeoMechanicsSoilStructureInteractionTests
from test_water_pressure import KratosGeoMechanicsWaterPressureTests
from test_dynamics import KratosGeoMechanicsDynamicsTests
from test_elements import KratosGeoMechanicsElementTypeTests
from test_steady_state_groundwater_flow import KratosGeoMechanicsSteadyStateGroundWaterFlowTests
from test_soil_weight import KratosGeoMechanicsSoilWeightTests


def AssambleTestSuites(is_team_city):
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    # Create an array with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample

    small_test_cases = [
        KratosGeoMechanicsGeneralTests,
        KratosGeoMechanicsExcavationTests,
        KratosGeoMechanicsInterfaceTests,
        KratosGeoMechanicsResetDisplacementTests,
        KratosGeoMechanicsSoilStructureInteractionTests,
        KratosGeoMechanicsWaterPressureTests,
        KratosGeoMechanicsBenchmarkSet1,
        KratosGeoMechanicsBenchmarkSet2,
        KratosGeoMechanicsElementTypeTests,
        KratosGeoMechanicsSteadyStateGroundWaterFlowTests,
        KratosGeoMechanicsSoilWeightTests
        ]

    # Create an array with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample

    night_test_cases = [KratosGeoMechanicsDynamicsTests]
    night_test_cases.extend(small_test_cases)

    # Create an array that contains all the tests from every testCase
    # in the list:

    all_test_cases = []
    all_test_cases.extend(night_test_cases)

    # add the tests to the corresponding suite,
    if is_team_city:

        smallSuite = unittest.TestSuite()
        nightSuite = unittest.TestSuite()
        allSuite = unittest.TestSuite()

        for test in small_test_cases:
            smallSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        for test in night_test_cases:
            nightSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        for test in all_test_cases:
            allSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        suites = allSuite
    else:
        suites = KratosUnittest.KratosSuites
        smallSuite = suites['small']
        nightSuite = suites['nightly']
        allSuite = suites['all']

        smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(small_test_cases))
        nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(night_test_cases))
        allSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(all_test_cases))

    return suites


if __name__ == '__main__':
    is_team_city = False

    try:
        from teamcity import is_running_under_teamcity
        from teamcity.unittestpy import TeamcityTestRunner

        is_team_city = is_running_under_teamcity()
    except ImportError:
        pass

    if is_team_city:
        import unittest
        runner = TeamcityTestRunner()
        runner.run(AssambleTestSuites(is_team_city))
    else:
        KratosUnittest.runTests(AssambleTestSuites(is_team_city))

