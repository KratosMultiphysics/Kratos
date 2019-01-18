from KratosMultiphysics import *
from KratosMultiphysics.MeshMovingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("MeshMovingApplicationFastSuite")

if __name__ == '__main__':
    run()
