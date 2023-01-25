import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication
from KratosMultiphysics import mpi

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    KratosMultiphysics.Tester.RunTestSuite("TrilinosApplicationFastSuite")

if __name__ == '__main__':
    run()