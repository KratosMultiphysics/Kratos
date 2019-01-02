from KratosMultiphysics import *
from KratosMultiphysics.ShallowWaterApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
    Tester.RunTestSuite("ShallowWaterApplicationFastSuite")

if __name__ == '__main__':
    run()
