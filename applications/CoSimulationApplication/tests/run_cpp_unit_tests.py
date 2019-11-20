from KratosMultiphysics import *
from KratosMultiphysics.CoSimulationApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosCoSimulationFastSuite")

if __name__ == '__main__':
    run()
