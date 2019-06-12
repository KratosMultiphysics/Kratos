from KratosMultiphysics import *
from KratosMultiphysics.ExaquteSandboxApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("ExaquteSandboxApplicationFastSuite")

if __name__ == '__main__':
    run()