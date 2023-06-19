# Definition of the classes for the SMALL TESTS
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ThermalDEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import TestFactory as TF

class SphereSphere_Touch_1(TF.TestFactory):
    file_name       = "Test_SphereSphere_Touch/SphereSphere_Touch"
    file_parameters = "Test_SphereSphere_Touch/ProjectParameters_1.json"

class SphereSphere_Touch_2(TF.TestFactory):
    file_name       = "Test_SphereSphere_Touch/SphereSphere_Touch"
    file_parameters = "Test_SphereSphere_Touch/ProjectParameters_2.json"

available_tests = [SphereSphere_Touch_1,SphereSphere_Touch_2]

def SetTestSuite(suites):
    small_suite = suites['small']
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))
    return small_suite
