from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    try:
        Tester.RunTestSuite("KratosSolidMechanicsFastSuite")
    except RuntimeError:
        print(" cpp test not included ")

if __name__ == '__main__':
    run()
