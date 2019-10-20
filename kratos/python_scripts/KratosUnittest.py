from __future__ import print_function, absolute_import, division
from KratosMultiphysics import Logger

from unittest import * # needed to make all functions available to the tests using this file
from unittest.util import safe_repr
from contextlib import contextmanager

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

    @classmethod
    def setUpClass(cls):
        if (sys.version_info < (3, 2)):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

    def run(self, result=None):
        super(TestCase,self).run(result)

    def failUnlessEqualWithTolerance(self, first, second, tolerance, msg=None):
        ''' fails if first and second have a difference greater than
        tolerance '''

        if first < (second - tolerance) or first > (second + tolerance):
            raise self.failureException(msg or '%r != %r within %r places' % (first, second, tolerance))

    def assertIsClose(self, first, second, rel_tol=None, abs_tol=None, msg=None):
        """Fail if the two objects are unequal as determined by their
           absolute and relative difference

           If the two objects compare equal then they will automatically
           compare relative almost equal.
        """
        if first == second:
            # shortcut
            return

        if rel_tol is None:
            rel_tol = 1e-09
        if abs_tol is None:
            abs_tol = 0.0

        if isclose(first, second, rel_tol, abs_tol):
            return

        standardMsg = '%s != %s within %s rel-tol and %s abs-tol' % (safe_repr(first),
                                                     safe_repr(second),
                                                     rel_tol, abs_tol)
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)


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
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\' \'(validation)\'',  # noqa
        '\t -v, --verbose: Verbosity level: 0, 1 (Default), 2',
        '\t --using-mpi: If running in MPI and executing the MPI-tests'
    ]
    for l in lines:
        Logger.PrintInfo(l) # using the logger to only print once in MPI

def main():
    # this deliberately overiddes the function "unittest.main",
    # because it cannot parse extra command line arguments
    if "--using-mpi" in sys.argv:
        sys.argv.remove("--using-mpi") # has to be removed bcs unittest cannot parse it
    import unittest
    unittest.main()

def runTests(tests):
    verbose_values = [0, 1, 2]
    level_values = ['all', 'small', 'nightly', 'validation']

    verbosity = 1
    level = 'all'
    is_mpi = False

    # Parse Commandline
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            'hv:l:', [
                'help',
                'verbose=',
                'level=',
                'using-mpi'
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
        elif o in ('--using-mpi'):
            is_mpi = True
        else:
            assert False, 'unhandled option'

    if is_mpi:
        level = "mpi_" + level

    if tests[level].countTestCases() == 0:
        print(
            '[Warning]: "{}" test suite is empty'.format(level),
            file=sys.stderr)
    else:
        result = not TextTestRunner(verbosity=verbosity, buffer=True).run(tests[level]).wasSuccessful()
        sys.exit(result)


KratosSuites = {
    'small':          TestSuite(),
    'nightly':        TestSuite(),
    'all':            TestSuite(),
    'validation':     TestSuite(),
    'mpi_small':      TestSuite(),
    'mpi_nightly':    TestSuite(),
    'mpi_all':        TestSuite(),
    'mpi_validation': TestSuite(),
}


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    '''Same implementation as math.isclose
    self-implemented bcs math.isclose was only introduced in python3.5
    see https://www.python.org/dev/peps/pep-0485/
    '''
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class WorkFolderScope:
    """ Helper-class to execute test in a specific target path
        Input
        -----
        - rel_path_work_folder: String
            Relative path of the target dir from the calling script

        - file_path: String
            Absolute path of the calling script

        - add_to_path: Bool
            "False" (default) if no need to add the target dir to the path, "True" otherwise.
    """

    def __init__(self, rel_path_work_folder, file_path, add_to_path=False):
        self.currentPath = os.getcwd()
        self.add_to_path = add_to_path
        if self.add_to_path:
            self.currentPythonpath = sys.path
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(file_path)), rel_path_work_folder))

    def __enter__(self):
        os.chdir(self.scope)
        if self.add_to_path:
            sys.path.append(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)
        if self.add_to_path:
            sys.path = self.currentPythonpath