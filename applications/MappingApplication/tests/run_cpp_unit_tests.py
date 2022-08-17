from KratosMultiphysics import Tester
import KratosMultiphysics.MappingApplication # registering the tests
import sys

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    if len(sys.argv) < 2:
        Tester.RunTestSuite("KratosMappingApplicationSerialTestSuite")
    else:
        Tester.RunTestCases(sys.argv[1])

if __name__ == '__main__':
    run()
