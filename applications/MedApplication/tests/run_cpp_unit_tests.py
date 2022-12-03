from KratosMultiphysics import *
from KratosMultiphysics.MedApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosMedFastSuite")

if __name__ == '__main__':
    run()
