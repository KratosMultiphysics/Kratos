import KratosMultiphysics as Kratos
import sys
def run():
    if len(sys.argv) < 2:
        Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)
        Kratos.Tester.RunAllTestCases()
    else :
        Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)
        Kratos.Tester.RunTestCases(sys.argv[1])

if __name__ == '__main__':
    run()