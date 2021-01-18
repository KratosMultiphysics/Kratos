from KratosMultiphysics import *
from KratosMultiphysics import mpi
from KratosMultiphysics.MappingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")

if __name__ == '__main__':
    run()
