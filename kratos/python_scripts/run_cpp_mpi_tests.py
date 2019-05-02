from KratosMultiphysics import *
from KratosMultiphysics import mpi

#Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)

Tester.RunTestSuite("KratosMPICoreFastSuite")