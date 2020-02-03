import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.CSharpWrapperApplication

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.TESTS_OUTPUTS) # TESTS_OUTPUTS
    KratosMultiphysics.Tester.RunTestSuite("KratosCSharpWrapperApplicationFastSuite")

if __name__ == '__main__':
    run()
