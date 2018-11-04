from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ContactStructuralMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosContactStructuralMechanicsFastSuite")

if __name__ == '__main__':
    run()
