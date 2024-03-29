from KratosMultiphysics import *
from KratosMultiphysics.MPMApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosMPMFastSuite")

if __name__ == '__main__':
    run()
