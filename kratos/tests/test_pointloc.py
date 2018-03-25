from KratosMultiphysics import *
# If you want to run tests defined in an application, import it here.

Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level

Tester.RunTestCases("*PointLocator*") #Test a specific case

# from KratosMultiphysics import *
# Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST)
# Tester.RunTestCases("*Serializer*")
