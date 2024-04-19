import os
import sys
import argparse

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.testing import utilities as testing_utils

def main():
    # Define the command
    cmd = testing_utils.GetPython3Command()

    # List of application
    applications = kratos_utils.GetListOfAvailableApplications()

    # Keep the worst exit code
    exit_code = 0

    # Parse Commandline
    # parse command line options
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--command', default=cmd, help="Use the provided command to launch test cases. If not provided, the default \'python\' executable is used")
    parser.add_argument('-l', '--level', default='all', choices=['all', 'nightly', 'small', 'validation'], help="Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'")
    parser.add_argument('-v', '--verbosity', default=1, type=int, choices=[0, 1, 2], help="Verbosity level: 0, 1 (Default), 2")
    parser.add_argument('-a', '--applications', default=applications, help="List of applications to run separated by \':\'. All compiled applications will be run by default")
    parser.add_argument('-m', '--using-mpi', default=False, help="If running in MPI and executing the MPI-tests")
    parser.add_argument('-t', '--timer', default=-1, help="Use the provided custom time limit for the execution. If not provided, the default values are used")

    args = parser.parse_args()

    # Set if mpi
    if args.using_mpi:
        level = "mpi_" + level

    # Parser the applications
    if isinstance(args.applications,str):
        parsedApps = args.applications.split(':')
    else:
        parsedApps = args.applications
    for a in parsedApps:
        if a not in applications + ['KratosCore']:
            print('Warning: Application {} does not exist'.format(a))
            sys.exit()
    applications = parsedApps
    if 'KratosCore' in applications:
        applications.remove('KratosCore')

    # Set timeout of the different levels
    signalTime = None
    if int(args.timer) > 0:
        signalTime = int(args.timer)
    else:
        if args.level == 'small':
            signalTime = int(90)
        elif args.level == 'nightly':
            signalTime = int(900)

    # Create the commands
    commander = testing_utils.Commander()

    exit_codes = {}

    # KratosCore must always be runned
    testing_utils.PrintTestHeader("KratosCore")

    with KratosUnittest.SupressConsoleOutput():
        commander.RunTestSuit(
            'KratosCore',
            'kratos',
            os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()),
            args.level,
            args.verbosity,
            cmd,
            signalTime
        )

    testing_utils.PrintTestFooter("KratosCore", commander.exitCode)
    exit_codes["KratosCore"] = commander.exitCode

    # Run the tests for the rest of the Applications
    for application in applications:
        testing_utils.PrintTestHeader(application)

        with KratosUnittest.SupressConsoleOutput():
            commander.RunTestSuit(
                application,
                application,
                KM.KratosPaths.kratos_applications+'/',
                args.level,
                args.verbosity,
                cmd,
                signalTime
            )

        testing_utils.PrintTestFooter(application, commander.exitCode)
        exit_codes[application] = commander.exitCode

    testing_utils.PrintTestSummary(exit_codes)
    sys.exit(max(exit_codes.values()))

if __name__ == "__main__":
    main()
