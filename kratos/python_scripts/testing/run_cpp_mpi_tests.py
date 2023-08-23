import argparse
import sys

from KratosMultiphysics import mpi, Tester

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--using-mpi',
        action='store_true',
        help="Force MPI mode if auto-detection fails")
    parser.add_argument(
        'match_string',
        metavar='MATCH-STRING',
        nargs='?',
        help="run cases with names matching MATCH-STRING")

    args = parser.parse_args()

    if args.match_string:
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
        return Tester.RunTestCases(args.match_string)
    else:
        Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
        return Tester.RunAllDistributedTestCases()

if __name__ == '__main__':
    try:
        exit_code = main()
    except Exception as e:
        print('[Error]:', e, file=sys.stderr)
        exit_code = 1

    sys.exit(exit_code)