from __future__ import print_function, absolute_import, division
from unittest import *

import getopt
import sys


class TestCase(TestCase):
    def failUnlessEqualWithTolerance(self, first, second, tolerance):
        if first < second + tolerance and first > second - tolerance:
            return True
        return False

    assertEqualTolerance = failUnlessEqualWithTolerance


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


def runTests(tests):
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

    TextTestRunner(verbosity=verbosity).run(tests[level])
