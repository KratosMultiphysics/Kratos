from __future__ import print_function, absolute_import, division

import imp
import os
import re
import getopt
import sys
import subprocess
import signal

from KratosMultiphysics import KratosLoader
from KratosMultiphysics.KratosUnittest import CaptureStdout


def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v verbosity] [-a app1:[app2:...]]',  # noqa
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t -a, --applications: List of applications to run separated by \':\'. All compiled applications will be run by default',  # noqa
        '\t -v, --verbose: Verbosity level: 0, 1 (Default), 2'
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


def RunTestSuit(application, applicationPath, path, level, verbose, command):
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
        print(script)

        if possiblePaths[0]['Found'] != possiblePaths[0]['Expected']:
            print(
                '[Warning]: Application has been found in "{}" directory but it was expected in "{}". Please check the naming convention.'.format(
                    possiblePaths[0]['Found'],
                    possiblePaths[0]['Expected']),
                file=sys.stderr)

        if os.path.isfile(script):
            subprocess.call([
                command,
                script,
                '-l'+level,
                '-v'+str(verbose)
            ])
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


def main():

    # We need to fetch the command who called us to avoid problems with
    # python versions such as running python3 while default is python2
    command = sys.executable

    verbose_values = [0, 1, 2]
    level_values = ['all', 'nightly', 'small', 'validation']

    # Set default values
    applications = GetAvailableApplication()
    verbosity = 1
    level = 'all'

    # Parse Commandline
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            'ha:v:l:', [
                'help',
                'applications=',
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
        elif o in ('-a', '--applications'):
            try:
                parsedApps = a.split(':')
            except:
                print('Error: Cannot parse applications.')
                Usage()
                sys.exit()

            for a in parsedApps:
                if a not in applications + ['KratosCore']:
                    print('Warning: Application {} does not exists'.format(a))
                    sys.exit()

            applications = parsedApps

            if 'KratosCore' in applications:
                applications.remove('KratosCore')

        else:
            assert False, 'unhandled option'

    # Capture stdout from KratosUnittest
    CaptureStdout()

    # Set a Timer
    signal.signal(signal.SIGALRM, handler)

    timedLevels = ['small', 'nightly']

    if level == 'small':
        signalTime = int(60)
    elif level == 'nightly':
        signalTime = int(900)

    # Define the command
    cmd = os.path.dirname(GetModulePath('KratosMultiphysics'))+'/'+'runkratos'

    # KratosCore must always be runned
    print('Running tests for KratosCore', file=sys.stderr)
    sys.stderr.flush()

    if level in timedLevels:
        signal.alarm(signalTime)
    try:
        RunTestSuit(
            'KratosCore',
            'kratos',
            os.path.dirname(GetModulePath('KratosMultiphysics')),
            level,
            verbosity,
            cmd
        )
    except Exception as exc:
        print('\nABORT: Tests for KratosCore took to long. Process Killed.', file=sys.stderr)

    # Run the tests for the rest of the Applications
    for application in applications:
        print('Running tests for {}'.format(application), file=sys.stderr)
        sys.stderr.flush()

        if level in timedLevels:
            signal.alarm(signalTime)

        try:
            RunTestSuit(
                application,
                application,
                KratosLoader.kratos_applications+'/',
                level,
                verbosity,
                cmd
            )
        except Exception as exc:
            print('\nABORT: Tests for {} took to long. Process Killed.'.format(application), file=sys.stderr)


if __name__ == "__main__":
    main()
