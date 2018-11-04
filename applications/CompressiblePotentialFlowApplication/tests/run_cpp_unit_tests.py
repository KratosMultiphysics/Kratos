from KratosMultiphysics import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("CompressiblePotentialApplicationFastSuite")

if __name__ == '__main__':
    run()
