from KratosMultiphysics import *
from KratosMultiphysics.HDF5Application import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosHDF5TestSuite")

if __name__ == '__main__':
    run()
