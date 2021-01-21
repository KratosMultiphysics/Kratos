import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics import kratos_utilities as kratos_utils

from KratosMultiphysics.testing import utilities as testing_utils

import sys, os
import argparse
from pathlib import Path

if KM.IsDistributedRun():
    raise Exception("cannot be run with MPI!")


def Usage():
    ''' Prints the usage of the script '''

    lines = [
        'Usage:',
        '\t python kratos_run_tests [-l level] [-v verbosity] [-a app1:[app2:...]]',  # noqa
        'Options',
        '\t -h, --help: Shows this command',
        '\t -l, --level: Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'',  # noqa
        '\t -a, --applications: List of applications to run separated by \':\'. All compiled applications will be run by default',  # noqa
        '\t -v, --verbose: Verbosity level: 0, 1 (Default), 2',
        '\t -c, --command: Use the provided command to launch test cases. If not provided, the default \'runkratos\' executable is used',
    ]
    for l in lines:
        print(l) # using the logger to only print once in MPI





def main():
    # Set default values
    applications = kratos_utils.GetListOfAvailableApplications()

    # parse command line options
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--command', default=testing_utils.GetPython3Command())
    parser.add_argument('-l', '--level', default='all', choices=['all', 'nightly', 'small', 'validation'])
    parser.add_argument('-v', '--verbosity', default=1, type=int, choices=[0, 1, 2])
    parser.add_argument('-a', '--applications', default=applications, choices=applications)
    parser.add_argument('-p', '--processes', type=int, default=4)

    try:
        args = parser.parse_args()
    except:
        Usage()
        sys.exit(2)

    # Keep the worst exit code
    exit_code = 0

    # Set timeout of the different levels
    signalTime = None
    if args.level == 'small':
        signalTime = 60
    elif args.level == 'nightly':
        signalTime = 900

    # Create the commands
    commander = testing_utils.Commander()

    exit_codes = {}

    # KratosMPICore must always be executed
    with KratosUnittest.SupressConsoleOutput():
        for np in range(2, args.processes+1):
            testing_utils.PrintTestHeader("KratosMPICore; {} processes".format(np))

            commander.RunMPITestSuit(
                'kratos/mpi',
                Path(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath())),
                np,
                args.level,
                args.verbosity,
                args.command,
                signalTime
            )

            testing_utils.PrintTestFooter("KratosMPICore; {} processes".format(np), commander.exitCode)
            exit_codes["KratosMPICore_{}_procs".format(np)] = commander.exitCode

    # Run the tests for the rest of the Applications
    for application in applications:
        for np in range(2, args.processes+1):
            testing_utils.PrintTestHeader("{}; {} processes".format(application, np))
            testing_utils.PrintTestHeader(application)

            with KratosUnittest.SupressConsoleOutput():
                commander.RunMPITestSuit(
                    application,
                    Path(KM.KratosPaths.kratos_applications),
                    np,
                    args.level,
                    args.verbosity,
                    args.command,
                    signalTime
                )

            testing_utils.PrintTestFooter("{}; {} processes".format(application, np), commander.exitCode)
            exit_codes["{}_{}_procs".format(application, np)] = commander.exitCode


    testing_utils.PrintTestSummary(exit_codes)
    sys.exit(max(exit_codes.values()))


if __name__ == "__main__":
    main()
