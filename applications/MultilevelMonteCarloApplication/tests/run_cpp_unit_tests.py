from KratosMultiphysics import *
from KratosMultiphysics.MultilevelMonteCarloApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosMultilevelMonteCarloApplicationFastSuite")

if __name__ == '__main__':
    run()