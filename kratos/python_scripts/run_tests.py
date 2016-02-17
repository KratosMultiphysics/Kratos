from __future__ import print_function, absolute_import, division

import imp
import os
import re
import getopt
import sys
import subprocess

from KratosMultiphysics import KratosLoader
from KratosMultiphysics.KratosUnittest import CaptureStdout


def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v vervosity] [-a app1:[app2:...]]',  # noqa
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t -a, --applications: List of applications to run separated by \':\'. All compiled applications will be run by default',  # noqa
        '\t -v, --verbose: Vervosty level: 0, 1 (Default), 2'
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
        f.split('.')[0] for f in os.listdir(kratosPath)
        if re.match(r'.*Application\.py$', f)
    ]

    return apps


def RunTestSuit(application, path, level, verbose):
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

    '''

    script = path+'/tests/'+'test_'+application+'.py'

    if os.path.isfile(script):
        subprocess.call([
            'python',
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
        if verbose > 1:
            print(
                '  expected file: "{}"'.format(
                    script),
                file=sys.stderr)


def main():

    verbose_values = [0, 1, 2]
    level_values = ['all', 'nightly', 'small']

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

    # KratosCore must always be runned
    print('Running tests for KratosCore')
    RunTestSuit(
        'KratosCore',
        os.path.dirname(GetModulePath('KratosMultiphysics'))+'/'+'kratos',
        level,
        verbosity
    )

    # Run the tests for the rest of the Applications
    for application in applications:
        print('Running tests for {}'.format(application))
        RunTestSuit(
            application,
            KratosLoader.kratos_applications+'/'+application,
            level,
            verbosity
        )

if __name__ == "__main__":
    main()
