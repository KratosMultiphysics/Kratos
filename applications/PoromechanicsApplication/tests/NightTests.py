# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class arc_length_test(TF.TestFactory):
     file_name = "element_tests/arc_length_test/arc_length_test"
     file_parameters = "element_tests/arc_length_test/ProjectParameters.json"
class consolidation_interface_3D(TF.TestFactory):
     file_name = "element_tests/consolidation_interface_3D/consolidation_interface_3D"
     file_parameters = "element_tests/consolidation_interface_3D/ProjectParameters.json"


def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            arc_length_test,
            consolidation_interface_3D
        ])
    )

    return night_suite
