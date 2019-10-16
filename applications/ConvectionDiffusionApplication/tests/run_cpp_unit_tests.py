from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosConvectionDiffusionFastSuite")

if __name__ == '__main__':
    run()
