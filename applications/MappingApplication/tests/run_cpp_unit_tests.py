from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosMappingApplicationSerialTestSuite")

if __name__ == '__main__':
    run()
