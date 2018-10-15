from __future__ import print_function, absolute_import, division

import imp
import os
import re
import getopt
import sys
import subprocess
import threading

from KratosMultiphysics import Tester
from KratosMultiphysics import KratosLoader
from KratosMultiphysics.KratosUnittest import SupressConsoleOutput, SupressConsoleError, SupressAllConsole


def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v verbosity] [-a app1:[app2:...]]',  # noqa
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t              For MPI tests, use the equivalent distributed test suites: \'(mpi_all)\', \'(mpi_nightly)\' \'(mpi_small)\'',
        '\t -a, --applications: List of applications to run separated by \':\'. All compiled applications will be run by default',  # noqa
        '\t -v, --verbose: Verbosity level: 0, 1 (Default), 2',
        '\t -c, --command: Use the provided command to launch test cases. If not provided, the default \'runkratos\' executable is used'
    ]

    for l in lines:
        print(l)


def GetModulePath(module):
    ''' Returns the location of a module using its absolute path

    Return
    ------
    string
        The absolute path of the module

    '''

    return imp.find_module(module)[1]


def GetAvailableApplication():
    ''' Return the list of applications available in KratosMultiphysics

    Return a list of compiled applications available in the KratosMultiphysics
    module.

    Return
    ------
    list of string
        List of the names of the applications

    '''
    kratosPath = GetModulePath('KratosMultiphysics')

    apps = [
        f.split('.')[0] for f in os.listdir(kratosPath) if re.match(r'.*Application\.py$', f)
    ]

    return apps


def handler(signum, frame):
    raise Exception("End of time")


class Commander(object):
    def __init__(self):
        self.process = None
        self.exitCode = 0

    def RunTestSuitInTime(self, application, applicationPath, path, level, verbose, command, timeout):
        if(timeout > -1):
            t = threading.Thread(
                target=self.RunTestSuit,
                args=(application, applicationPath, path, level, verbose, command)
            )

            t.start()
            t.join(timeout)

            if t.isAlive():
                self.process.terminate()
                t.join()
                print('\nABORT: Tests for {} took to long. Process Killed.'.format(application), file=sys.stderr)
            else:
                print('\nTests for {} finished in time ({}s).'.format(application, timeout))
        else:
            self.RunTestSuit(application, applicationPath, path, level, verbose, command)

    def RunTestSuit(self, application, applicationPath, path, level, verbose, command):
        ''' Calls the script that will run the tests.

        Input
        -----
        application: string
            Name of the application that will be tested.

        path: string
            Absoulte path with the location of the application.

        level: string
            minimum level of the test that will be run if possible.

        verbose: int
            detail of the ouptut. The grater the verbosity level, the greate the
            detail will be.

        command: string
            command to be used to call the tests. Ex: Python, Python3, Runkratos

        '''

        self.exitCode = 0
        appNormalizedPath = applicationPath.lower().replace('_', '')

        possiblePaths = [
            {'Found': p, 'FoundNormalized': p.split('/')[-1].lower().replace('_', ''), 'Expected': applicationPath, 'ExpectedNormalized': appNormalizedPath} for p in os.listdir(path) if p.split('/')[-1].lower().replace('_', '') == appNormalizedPath
        ]

        if len(possiblePaths) < 1:
            if verbose > 0:
                print(
                    '[Warning]: No directory found for {}'.format(
                        application),
                    file=sys.stderr)
                sys.stderr.flush()
        elif len(possiblePaths) > 1:
            if verbose > 0:
                print('Unable to determine correct path for {}'.format(application), file=sys.stderr)
                print(
                    'Please try to follow the standard naming convention \'FooApplication\' Snake-Capital string  without symbols.',
                    file=sys.stderr)
            if verbose > 1:
                print('Several possible options were found:', file=sys.stderr)
                for p in possiblePaths:
                    print('\t', p, file=sys.stderr)
        else:
            script = path+'/'+possiblePaths[0]['Found']+'/tests/'+'test_'+application+'.py'
            print(script, file=sys.stderr)

            if possiblePaths[0]['Found'] != possiblePaths[0]['Expected']:
                print(
                    '[Warning]: Application has been found in "{}" directory but it was expected in "{}". Please check the naming convention.'.format(
                        possiblePaths[0]['Found'],
                        possiblePaths[0]['Expected']),
                    file=sys.stderr)

            if os.path.isfile(script):
                self.process = subprocess.Popen([
                    command,
                    script,
                    '-l'+level,
                    '-v'+str(verbose)
                ], stdout=subprocess.PIPE)

                # Used instead of wait to "soft-block" the process and prevent deadlocks
                # and capture the first exit code different from OK
                process_stdout, process_stderr = self.process.communicate()
                if process_stdout:
                    print(process_stdout.decode('ascii'), file=sys.stdout)
                if process_stderr:
                    print(process_stderr.decode('ascii'), file=sys.stderr)

                # Running out of time in the tests will send the error code -15. We may want to skip
                # that one in a future. Right now will throw everything different from 0.
                self.exitCode = int(self.process.returncode != 0)
            else:
                if verbose > 0:
                    print(
                        '[Warning]: No test script found for {}'.format(
                            application),
                        file=sys.stderr)
                    sys.stderr.flush()
                if verbose > 1:
                    print(
                        '  expected file: "{}"'.format(
                            script),
                        file=sys.stderr)
                    sys.stderr.flush()

    def RunCppTests(self, applications):
        ''' Calls the cpp tests directly
        '''

        self.exitCode = 0

        # importing the apps such that they get registered for the cpp-tests
        for application in applications:
            __import__("KratosMultiphysics." + application)

        try:
            Tester.SetVerbosity(Tester.Verbosity.PROGRESS)
            self.exitCode = Tester.RunAllTestCases()
        except Exception as e:
            print('[Warning]:', e, file=sys.stderr)
            self.exitCode = 1


def main():

    # Define the command
    cmd = os.path.dirname(GetModulePath('KratosMultiphysics'))+'/'+'runkratos'

    verbose_values = [0, 1, 2]
    level_values = ['all', 'nightly', 'small', 'validation', 'mpi_all', 'mpi_small', 'mpi_nightly', 'mpi_validation']

    # Set default values
    applications = GetAvailableApplication()
    verbosity = 1
    level = 'all'

    # Keep the worst exit code
    exit_code = 0

    # Parse Commandline
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            'ha:v:l:c:', [
                'help',
                'applications=',
                'verbose=',
                'level=',
                'command='
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
        elif o in ('-a', '--applications'):
            try:
                parsedApps = a.split(':')
            except:
                print('Error: Cannot parse applications.')
                Usage()
                sys.exit()

            for a in parsedApps:
                if a not in applications + ['KratosCore']:
                    print('Warning: Application {} does not exist'.format(a))
                    sys.exit()

            applications = parsedApps

            if 'KratosCore' in applications:
                applications.remove('KratosCore')

        elif o in ('-c', '--command'):
            try:
                cmd = str(a)
            except:
                print('Error: Cannot parse command name {0}.'.format(a))
                Usage()
                sys.exit()

        else:
            assert False, 'unhandled option'

    # Set timeout of the different levels
    signalTime = int(-1)
    if level == 'small':
        signalTime = int(60)
    elif level == 'nightly':
        signalTime = int(900)

    # Create the commands
    commander = Commander()

    # KratosCore must always be runned
    print('Running tests for KratosCore', file=sys.stderr)

    with SupressConsoleOutput():
        commander.RunTestSuitInTime(
            'KratosCore',
            'kratos',
            os.path.dirname(GetModulePath('KratosMultiphysics')),
            level,
            verbosity,
            cmd,
            signalTime
        )

    exit_code = max(exit_code, commander.exitCode)

    # Run the tests for the rest of the Applications
    for application in applications:
        print('Running tests for {}'.format(application), file=sys.stderr)
        sys.stderr.flush()

        with SupressConsoleOutput():
            commander.RunTestSuitInTime(
                application,
                application,
                KratosLoader.kratos_applications+'/',
                level,
                verbosity,
                cmd,
                signalTime
            )

        exit_code = max(exit_code, commander.exitCode)

    # Run the cpp tests (does the same as run_cpp_tests.py)
    print('Running cpp tests', file=sys.stderr)
    with SupressConsoleOutput():
        commander.RunCppTests(applications)

    exit_code = max(exit_code, commander.exitCode)

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
