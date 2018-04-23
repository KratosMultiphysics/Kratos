from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosStructuralMechanicsFastSuite")

if __name__ == '__main__':
    run()
