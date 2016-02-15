from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_model_part_io import TestModelPartIO as TModelPartIO
from test_model_part import TestModelPart as TModelPart


def AssambleTestSuits():
    ''' Generates the test suits to run.

    Generates the test suits to run. At least, it needs to return the suits:
    "small", "nighlty" and "all"

    Return
    ------

    set of suits:
        A set of suits in the form of {"suitname": suit}.

    '''

    # Create a test suit with the selected tests (Small tests):
    SmallSuite = KratosUnittest.TestSuite()
    SmallSuite.addTest(TModelPartIO('test_model_part_io_read_model_part'))
    SmallSuite.addTest(TModelPart('test_model_part_properties'))

    # Create a test suit with the selected tests plus all small tests
    NightSuite = KratosUnittest.TestSuite([SmallSuite])
    NightSuite.addTest(TModelPart('test_model_part_sub_model_parts'))
    NightSuite.addTest(TModelPart('test_model_part_nodes'))
    NightSuite.addTest(TModelPart('test_model_part_tables'))

    # Create a test suit that contains all the tests:
    AllModelPartIO = KratosUnittest.TestLoader().loadTestsFromTestCase(
        TModelPartIO)
    AllModelPart = KratosUnittest.TestLoader().loadTestsFromTestCase(
        TModelPart)

    AllSuit = KratosUnittest.TestSuite([
        AllModelPartIO,
        AllModelPart
    ])

    return {
        'small': SmallSuite,
        'nightly': NightSuite,
        'all': AllSuit
    }

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuits())
