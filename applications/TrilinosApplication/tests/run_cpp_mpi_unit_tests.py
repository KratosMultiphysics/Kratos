import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication
import sys

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    if len(sys.argv) < 2:
        KratosMultiphysics.Tester.RunTestSuite("KratosTrilinosApplicationMPITestSuite")
    else:
        KratosMultiphysics.Tester.RunTestCases(sys.argv[1])

if __name__ == '__main__':
    run()