from KratosMultiphysics import *
from KratosMultiphysics.ExaquteSandboxApplication import *
try:
    from KratosMultiphysics.MeshingApplication import *
except:
    raise Exception ("MeshingApplication is not compiled but is required to run the ExaquteSandboxApplication tests")

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("ExaquteSandboxApplicationFastSuite")

if __name__ == '__main__':
    run()
