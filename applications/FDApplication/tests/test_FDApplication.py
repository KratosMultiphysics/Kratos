# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FDApplication import *

import sys

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from SmallTests import Bfecc as TBfecc


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
    SmallSuite.addTest(TBfecc('testCreateSolver'))

    # Create a test suit with the selected tests plus all small tests
    NightSuite = KratosUnittest.TestSuite([SmallSuite])

    # Create a test suit that contains all the tests:
    AllBfecc = KratosUnittest.TestLoader().loadTestsFromTestCase(
        TBfecc)

    AllSuit = KratosUnittest.TestSuite([
        AllBfecc,
    ])

    return {
        'small': SmallSuite,
        'nightly': NightSuite,
        'all': AllSuit
    }

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuits())
