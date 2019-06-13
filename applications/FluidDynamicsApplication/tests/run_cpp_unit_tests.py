from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("FluidDynamicsApplicationFastSuite")

if __name__ == '__main__':
    run()
