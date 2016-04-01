from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites
from test_kratos_parameters import TestParameters as TParameters
from test_model_part_io import TestModelPartIO as TModelPartIO
from test_model_part import TestModelPart as TModelPart


def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']

    smallSuite.addTest(TModelPartIO('test_model_part_io_read_model_part'))
    smallSuite.addTest(TModelPart('test_model_part_properties'))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']

    nightSuite.addTests(map(TModelPart, [
        'test_model_part_sub_model_parts',
        'test_model_part_nodes',
        'test_model_part_tables'
    ]))

    nightSuite.addTests(map(TParameters, [
        'test_kratos_parameters',
        'test_kratos_change_parameters',
        'test_kratos_copy_parameters',
        'test_kratos_wrong_parameters'
    ]))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TModelPartIO,
            TModelPart,
            TParameters
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
