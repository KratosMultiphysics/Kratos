from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosMappingApplicationTestSuite")

if __name__ == '__main__':
    run()
