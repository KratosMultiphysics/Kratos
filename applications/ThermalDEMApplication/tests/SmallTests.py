# Definition of the classes for the SMALL TESTS
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ThermalDEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import TestFactory as TF

class test(TF.TestFactory):
    file_name       = "test"
    file_parameters = "ProjectParameters.json"

available_tests = [test]

def SetTestSuite(suites):
    small_suite = suites['small']
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))
    return small_suite
