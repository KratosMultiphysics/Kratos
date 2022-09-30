import KratosMultiphysics as Kratos
import sys
if len(sys.argv) < 2:
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_LIST)
    Kratos.Tester.RunAllTestCases()
else :
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)
    Kratos.Tester.RunTestCases(sys.argv[1])

