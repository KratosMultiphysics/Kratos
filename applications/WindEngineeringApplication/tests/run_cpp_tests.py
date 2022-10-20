from KratosMultiphysics import *
from KratosMultiphysics.WindEngineeringApplication import *

def Run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS)
    Tester.RunTestSuite("KratosWindEngineeringFastSuite")

if __name__ == '__main__':
    Run()
