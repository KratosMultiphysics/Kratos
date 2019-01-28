from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    try:
        Tester.RunTestSuite("KratosSolidMechanicsFastSuite")
    except:
        print(" cpp test not included ")
        pass

if __name__ == '__main__':
    run()
