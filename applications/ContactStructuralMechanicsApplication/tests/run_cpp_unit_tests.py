import KratosMultiphysics
import KratosMultiphysics.ContactStructuralMechanicsApplication

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    KratosMultiphysics.Tester.RunTestSuite("KratosContactStructuralMechanicsFastSuite")

if __name__ == '__main__':
    run()
