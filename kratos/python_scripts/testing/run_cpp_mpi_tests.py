import sys

from KratosMultiphysics import mpi, Tester

if __name__ == '__main__':
    if len(sys.argv) < 2:
        Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
        Tester.RunAllDistributedTestCases()
    else :
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
        Tester.RunTestCases(sys.argv[1])
