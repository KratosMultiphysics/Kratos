# import Kratos

import os
import sys
import argparse

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
from test_transient_groundwater_flow import KratosGeoMechanicsTransientGroundWaterFlowTests
from test_soil_weight import KratosGeoMechanicsSoilWeightTests
from test_settlement import KratosGeoMechanicsSettlementTests
from test_curved_beam_elements import KratosGeoMechanicsCurvedBeamElementTests
from test_elementary_groundwater_flow import TestElementaryGroundWaterFlow
from test_sellmeijers_rule import TestSellmeijersRule
from test_consecutive_pipe_lines import TestConsecutivePipeLines

# cpp tests
from test_piping_element_unit import TestUnitPipingElements
from test_normal_flux_condition import TestNormalFluxCondition


def AssembleTestSuites(is_team_city):
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

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
        KratosGeoMechanicsTransientGroundWaterFlowTests,
        KratosGeoMechanicsSoilWeightTests,
        KratosGeoMechanicsSettlementTests,
        KratosGeoMechanicsCurvedBeamElementTests,
        TestNormalFluxCondition
    ]

    # Create an array with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample

    night_test_cases = [KratosGeoMechanicsDynamicsTests,
                        TestSellmeijersRule,
                        TestElementaryGroundWaterFlow]

    # Create an array with all long tests only for validations
    
    valid_test_cases = [TestConsecutivePipeLines]
    
    
    # Create an array that contains all the tests from every testCase
    # in the list:

    all_test_cases = []
    all_test_cases.extend(night_test_cases)
    all_test_cases.extend(small_test_cases)
    all_test_cases.extend(valid_test_cases)
    suites = KratosUnittest.KratosSuites

    # add the tests to the corresponding suite,
    if is_team_city:
        smallSuite = unittest.TestSuite()
        nightSuite = unittest.TestSuite()
        validSuite = unittest.TestSuite()
        allSuite = unittest.TestSuite()

        small_test_cases.append(TestUnitPipingElements)

        for test in small_test_cases:
            smallSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        for test in night_test_cases:
            nightSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))
                
        
        for test in valid_test_cases:
            validSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        for test in all_test_cases:
            allSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(
                test))

        # suites = allSuite
        suites['small'] = smallSuite
        suites['nightly'] = nightSuite
        suites['validation'] = validSuite
        suites['all'] = allSuite
    else:
        smallSuite = suites['small']
        nightSuite = suites['nightly']
        validSuite = suites['validation']
        allSuite = suites['all']

        smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(small_test_cases))
        night_test_cases.extend(small_test_cases)
        nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(night_test_cases))
        validSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(valid_test_cases))
        allSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(all_test_cases))

    return suites


def get_level_argument():
    """
    This function is used only when the unit tests are run under teamcity.
    In this case the level the KratosUnittest class is not directly used to determine the level of unit tests
    Therefore we implement a simple function that enables this functionality

    :return: level: str
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--level', default='all', choices=['all', 'nightly', 'small', 'validation'])
    args = parser.parse_args()
    return args.level


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

        level = get_level_argument()
        test_suite_dictionary = AssembleTestSuites(is_team_city)
        tests_to_run = test_suite_dictionary[level]
        runner = TeamcityTestRunner()
        runner.run(tests_to_run)
    else:
        KratosUnittest.runTests(AssembleTestSuites(is_team_city))
