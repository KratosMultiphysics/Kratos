from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosSolidMechanicsFastSuite")

if __name__ == '__main__':
    run()
