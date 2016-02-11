from __future__ import print_function, absolute_import, division

import getopt
import sys

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_model_part_io import TestModelPartIO as TModelPartIO
from test_model_part import TestModelPart as TModelPart


def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v vervosity]',
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t -v, --verbose: Vervosty level: 0, 1 (Default), 2'
    ]

    for l in lines:
        print(l)


def AssambleTestSuits():
    ''' Generates the test suits to run.

    Generates the test suits to run. At least it will return the suits:
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

    tests = AssambleTestSuits()

    verbose_values = [0, 1, 2]
    level_values = ['all', 'small', 'nightly']

    verbosity = 1
    level = 'all'

    # Parse Commandline
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            'hv:l:', [
                'help',
                'verbose=',
                'level='
            ])
    except getopt.GetoptError as err:
        print(str(err))
        Usage()
        sys.exit(2)

    for o, a in opts:
        if o in ('-v', '--verbose'):
            if int(a) in verbose_values:
                verbosity = int(a)
            else:
                print('Error: {} is not a valid verbose level.'.format(a))
                Usage()
                sys.exit()
        elif o in ('-h', '--help'):
            Usage()
            sys.exit()
        elif o in ('-l', '--level'):
            if a in level_values:
                level = a
            else:
                print('Error: {} is not a valid level.'.format(a))
                Usage()
                sys.exit()
        else:
            assert False, 'unhandled option'

    KratosUnittest.TextTestRunner(verbosity=verbosity).run(tests[level])
