from KratosMultiphysics import *
from KratosMultiphysics.ExaquteSandboxApplication import *
try:
    from KratosMultiphysics.MeshingApplication import *
except ImportError:
    raise ImportError("Meshing application not imported. Must be compiled in order to run the tests.")

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("ExaquteSandboxApplicationFastSuite")

if __name__ == '__main__':
    run()