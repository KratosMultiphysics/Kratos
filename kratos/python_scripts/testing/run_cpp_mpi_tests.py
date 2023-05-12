import sys

from KratosMultiphysics import mpi, Tester

if __name__ == '__main__':
    #Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)

    if len(sys.argv) < 2:
        Tester.RunTestSuite("KratosMPICoreFastSuite")
    else:
        Tester.RunTestCases(sys.argv[1])
