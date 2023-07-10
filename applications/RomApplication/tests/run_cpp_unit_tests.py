from KratosMultiphysics import *
from KratosMultiphysics.RomApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS PROGRESS
    Tester.RunTestSuite("RomApplicationFastSuite")

if __name__ == '__main__':
    run()