from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosIgaFastSuite")

if __name__ == '__main__':
    run()
