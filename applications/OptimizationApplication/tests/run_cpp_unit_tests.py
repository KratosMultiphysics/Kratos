from KratosMultiphysics import *
from KratosMultiphysics.OptimizationApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosOptimizationFastSuite")

if __name__ == '__main__':
    run()
