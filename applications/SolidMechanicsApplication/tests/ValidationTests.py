# Definition of the classes for the VALIDATION TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest


#class ValidationTest(TF.TestFactory):
    #file_name = "path_to_my_test"


def SetTestSuite(suites):
    validation_suite = suites['validation']

    return validation_suite
