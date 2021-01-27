import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics import kratos_utilities as kratos_utils

from KratosMultiphysics.testing import utilities as testing_utils

import sys, os
import argparse
from pathlib import Path

if KM.IsDistributedRun():
    raise Exception("cannot be run with MPI!")


def main():
    # Set default values
    applications = kratos_utils.GetListOfAvailableApplications()

    # parse command line options
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--command', default=testing_utils.GetPython3Command())
    parser.add_argument('-l', '--level', default='all', choices=['all', 'nightly', 'small', 'validation'])
    parser.add_argument('-v', '--verbosity', default=1, type=int, choices=[0, 1, 2])
    # parser.add_argument('-a', '--applications', default=applications, choices=applications) # TODO
    parser.add_argument('-n', '--processes', type=int, default=4)
    parser.add_argument('-m', '--mpi_command', default="mpiexec")
    parser.add_argument('-f', '--mpi_flags', default="")
    parser.add_argument('-p', '--num_processes_flag', default="-np")

    args = parser.parse_args()

    # Set timeout of the different levels
    signalTime = None
    if args.level == 'small':
        signalTime = 60
    elif args.level == 'nightly':
        signalTime = 900

    # Create the commands
    commander = testing_utils.Commander()

    exit_codes = {}

    testing_utils.PrintTestHeader("KratosMPICore")
    # KratosMPICore must always be executed
    with KratosUnittest.SupressConsoleOutput():
        commander.RunMPITestSuit(
            'KratosMPICore',
            Path(os.path.dirname(kratos_utils.GetKratosMultiphysicsPath()))/"kratos"/"mpi",
            args.mpi_command,
            args.mpi_flags,
            args.num_processes_flag,
            args.processes,
            args.level,
            args.verbosity,
            args.command,
            signalTime
        )

    testing_utils.PrintTestFooter("KratosMPICore", commander.exitCode)
    exit_codes["KratosMPICore"] = commander.exitCode

    # Run the tests for the rest of the Applications
    for application in applications:
        testing_utils.PrintTestHeader(application)

        with KratosUnittest.SupressConsoleOutput():
            commander.RunMPITestSuit(
                application+"_mpi",
                Path(KM.KratosPaths.kratos_applications) / application,
                args.mpi_command,
                args.mpi_flags,
                args.num_processes_flag,
                args.processes,
                args.level,
                args.verbosity,
                args.command,
                signalTime
            )

        testing_utils.PrintTestFooter(application, commander.exitCode)
        exit_codes[application] = commander.exitCode

    testing_utils.PrintTestSummary(exit_codes)
    sys.exit(max(exit_codes.values()))


if __name__ == "__main__":
    main()
