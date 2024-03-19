import argparse
import sys
import importlib

from KratosMultiphysics import Tester

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'match_string',
        metavar='MATCH-STRING',
        nargs='?',
        help="run cases with names matching MATCH-STRING")
    parser.add_argument(
        '--run-all-apps',
        default=False,
        action='store_true',
        help="If all applicaions tests are run, not only core")
    
    args = parser.parse_args()

    # Imporing all available applications to make sure that the tests are registered
    if args.run_all_apps:
        from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications
        available_apps = GetListOfAvailableApplications()
        for app in available_apps:
            importlib.import_module("KratosMultiphysics." + app)

    if args.match_string:
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
        return Tester.RunTestCases(args.match_string)
    else:
        Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
        return Tester.RunAllTestCases()

if __name__ == '__main__':
    try:
        exit_code = main()
    except Exception as e:
        print('[Error]:', e, file=sys.stderr)
        exit_code = 1

    sys.exit(exit_code)