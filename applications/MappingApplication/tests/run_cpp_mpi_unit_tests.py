from KratosMultiphysics import Tester
from KratosMultiphysics.MappingApplication import MPIExtension # registering the tests
import sys

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    if len(sys.argv) < 2:
        Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")
    else:
        Tester.RunTestCases(sys.argv[1])

if __name__ == '__main__':
    run()
