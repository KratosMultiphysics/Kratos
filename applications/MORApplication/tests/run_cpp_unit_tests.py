from KratosMultiphysics import *
from KratosMultiphysics.MORApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # TESTS_OUTPUTS  PROGRESS
    Tester.RunTestSuite("KratosMORFastSuite")

if __name__ == '__main__':
    run()
