from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
try:
  from KratosMultiphysics.StructuralMechanicsApplication import *
except ImportError as e:
    Logger.PrintWarning("CPP Unittests", "StructuralMechanicsApplication not compiled")

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosMeshingApplicationFastSuite")

if __name__ == '__main__':
    run()
