import KratosMultiphysics as KM
from KratosMultiphysics.MedApplication import *
import sys

if len(sys.argv) < 2:
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    KM.Tester.RunTestSuite("KratosMedFastSuite")
else:
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.TESTS_OUTPUTS)
    KM.Tester.RunTestCases(sys.argv[1])
