import sys
sys.path.append(R"C:\software\Kratos")

from KratosMultiphysics import *
Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
Tester.RunAllTestCases()
