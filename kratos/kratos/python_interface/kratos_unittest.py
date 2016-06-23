from __future__ import print_function, absolute_import, division
from unittest import *

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

    def failUnlessEqualWithTolerance(self, first, second, tolerance):
        ''' fails if first and second have a difference greater than
        tolerance '''

        if first < second + tolerance and first > second - tolerance:
            return True
        return False

    assertEqualTolerance = failUnlessEqualWithTolerance


def CaptureStdout(bufferFile=None):
    ''' Captures stdout and redirects it to bufferFile. If no bufferFile
    is provided stdout is redirected to os.devnull by default '''

    sys.stdout.flush()
    newstdout = os.dup(1)

    if bufferFile is None:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, 1)
        os.close(devnull)
    else:
        os.dup2(bufferFile, 1)

    sys.stdout = os.fdopen(newstdout, 'w')


def CaptureStderr(bufferFile=None):
    ''' Captures stderr and redirects it to bufferFile. If no bufferFile
    is provided stderr is redirected to os.devnull by default '''

    sys.stderr.flush()
    newstderr = os.dup(2)

    if bufferFile is None:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, 2)
        os.close(devnull)
    else:
        os.dup2(bufferFile, 2)

    sys.stderr = os.fdopen(newstderr, 'w')


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
    level_values = ['all', 'small', 'nightly', 'validation']

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
        TextTestRunner(verbosity=verbosity, buffer=True).run(tests[level])

KratosSuites = {
    'small': TestSuite(),
    'nightly': TestSuite(),
    'all': TestSuite(),
    'validation': TestSuite()
}
