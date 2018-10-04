from __future__ import print_function, absolute_import, division
from unittest import *
from contextlib import contextmanager
import KratosMultiphysics

import getopt
import sys
import os


class TestLoader(TestLoader):
    def loadTestsFromTestCases(self, testCaseClasses):
        ''' Return a list of suites with all tests cases contained in every
        testCaseClass in testCaseClasses '''

        allTests = []

        for caseClasses in testCaseClasses:
            caseTests = self.loadTestsFromTestCase(caseClasses)
            allTests.append(caseTests)

        return allTests


class TestCase(TestCase):
    
    def run(self, result=None): 
        super(TestCase,self).run(result)        

    def failUnlessEqualWithTolerance(self, first, second, tolerance, msg=None):
        ''' fails if first and second have a difference greater than
        tolerance '''

        if first < (second - tolerance) or first > (second + tolerance):
            raise self.failureException(msg or '%r != %r within %r places' % (first, second, tolerance))

    assertEqualTolerance = failUnlessEqualWithTolerance

@contextmanager
def SupressConsoleOutput():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

@contextmanager
def SupressConsoleError():
    with open(os.devnull, "w") as devnull:
        old_stderr = sys.stderr
        sys.stderr = devnull
        try:  
            yield
        finally:
            sys.stderr = old_stderr

@contextmanager
def SupressAllConsole():
    with open(os.devnull, "w") as devnull:
        old_stderr = sys.stderr
        old_stdout = sys.stdout
        sys.stderr = devnull
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stderr = old_stderr
            sys.stdout = old_stdout

def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v verbosity]',
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t              For MPI tests, use the equivalent distributed test suites: \'(mpi_all)\', \'(mpi_nightly)\' \'(mpi_small)\'',
        '\t -v, --verbose: Verbosity level: 0, 1 (Default), 2'
    ]

    for l in lines:
        print(l)


def runTests(tests):
    verbose_values = [0, 1, 2]
    level_values = ['all', 'small', 'nightly', 'validation', 'mpi_all', 'mpi_small', 'mpi_nightly', 'mpi_validation']

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

    if tests[level].countTestCases() == 0:
        print(
            '[Warning]: "{}" test suite is empty'.format(level),
            file=sys.stderr)
    else:
        result = not TextTestRunner(verbosity=verbosity, buffer=True).run(tests[level]).wasSuccessful()
        sys.exit(result)


KratosSuites = {
    'small': TestSuite(),
    'nightly': TestSuite(),
    'all': TestSuite(),
    'validation': TestSuite(),
    'mpi_small': TestSuite(),
    'mpi_nightly': TestSuite(),
    'mpi_all': TestSuite(),
    'mpi_validation': TestSuite(),
}
