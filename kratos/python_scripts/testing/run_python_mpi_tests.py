import sys
import argparse
import multiprocessing

from pathlib import Path

import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities as kratos_utils

from KratosMultiphysics.testing import utilities as testing_utils

if KM.IsDistributedRun():
    raise Exception("cannot be run with MPI!")

def main():
    """
    This script executes the MPI tests.
    It takes no parameters and returns no values. The variables that are set by this function, as well as the tests
    that are executed, are determined by the command line arguments passed to the function. These arguments include:

    -c, --command: The command to use to launch test cases. If not provided, the default 'python' executable is used.
    -l, --level: The minimum level of detail of the tests. Choices are 'mpi_all' (default), 'mpi_nightly', 'mpi_small', and 'mpi_validation'.
    -v, --verbosity: The verbosity level. Choices are 0, 1 (default), and 2.
    -a, --applications: A list of applications to run, separated by ':'. All compiled applications will be run by default.
    -n, --processes: The number of processes to consider. Default is the number of cores of the system.
    -m, --mpi_command: The MPI command to consider. Default is mpiexec.
    -f, --mpi_flags: The additional MPI flags to consider. Default is empty.
    -p, --num_processes_flag: The flag used in order to introduce the number of processes considered.
    -t, --timer: The custom time limit for the execution. If not provided, the default values are used.

    In addition to setting default values and parsing command line arguments, the function also runs tests for 
    the 'KratosMPICore' application and for any other applications specified by the user. The tests are run using 
    the specified command and number of processes, and the results are printed to the console.
    """

    # Define the command
    cmd = testing_utils.GetPython3Command()

    # Set default values
    applications = kratos_utils.GetListOfAvailableApplications() + ["KratosMPICore"]

    # parse command line options
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--command', default=testing_utils.GetPython3Command(), help="Use the provided command to launch test cases. If not provided, the default \'python\' executable is used")
    parser.add_argument('-l', '--level', default='mpi_all', choices=['mpi_all', 'mpi_nightly', 'mpi_small', 'mpi_validation'], help="Minimum level of detail of the tests: \'all\'(Default) \'(nightly)\' \'(small)\'")
    parser.add_argument('-v', '--verbosity', default=1, type=int, choices=[0, 1, 2], help="Verbosity level: 0, 1 (Default), 2")
    parser.add_argument('-a', '--applications', default=applications, help="List of applications to run separated by \':\'. All compiled applications will be run by default")
    parser.add_argument('-n', '--processes', type=int, default=multiprocessing.cpu_count(), help="Number of processes considered. Default is the number of cores of the system")
    parser.add_argument('-m', '--mpi_command', default="mpiexec", help="MPI command considered. Default is mpiexec")
    parser.add_argument('-f', '--mpi_flags', default="", help="The additional MPI flags considered. Default is empty")
    parser.add_argument('-p', '--num_processes_flag', default="-np", help="Flag used in order to introduce the number of processes considered")
    parser.add_argument('-t', '--timer', default=-1, help="Use the provided custom time limit for the execution. If not provided, the default values are used")

    args = parser.parse_args()

    # Parser the applications
    if isinstance(args.applications,str):
        parsedApps = args.applications.split(':')
    else:
        parsedApps = args.applications
    for a in parsedApps:
        if a not in applications:
            print('Warning: Application {} does not exist'.format(a))
            sys.exit()
    applications = parsedApps

    # Set timeout of the different levels
    signalTime = None
    if int(args.timer) > 0:
        signalTime = int(args.timer)
    else:
        if args.level == 'mpi_small':
            signalTime = int(90)
        elif args.level == 'mpi_nightly':
            signalTime = int(900)

    # Create the commands
    commander = testing_utils.Commander()

    # Run the tests
    commander.RunMPIPythonTests(applications, args.mpi_command, args.mpi_flags, args.num_processes_flag, args.processes, args.level, args.verbosity, cmd, signalTime)

    # Exit message
    testing_utils.PrintTestSummary(commander.exitCodes)

    # Propagate exit code and end
    try:
        exit_code = max(commander.exitCodes.values())
    except:
        print("Failed to run tests")
        exit_code = 1

    sys.exit(exit_code)

if __name__ == "__main__":
    main()
