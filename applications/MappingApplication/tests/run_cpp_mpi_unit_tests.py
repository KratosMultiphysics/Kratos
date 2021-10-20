from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import MPIExtension # registering the tests

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")

if __name__ == '__main__':
    run()
