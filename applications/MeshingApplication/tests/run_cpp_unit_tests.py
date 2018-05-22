from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosMeshingApplicationFastSuite")

if __name__ == '__main__':
    run()
