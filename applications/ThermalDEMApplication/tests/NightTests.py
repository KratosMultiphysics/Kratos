# Definition of the classes for the NIGHTLY TESTS
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ThermalDEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import TestFactory as TF

available_tests = []

def SetTestSuite(suites):
    night_suite = suites['nightly']
    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))
    return night_suite
