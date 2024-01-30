# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class undrained_soil_column_2D(TF.TestFactory):
    file_name = "element_tests/undrained_soil_column_2D/undrained_soil_column_2D"
    file_parameters = "element_tests/undrained_soil_column_2D/ProjectParameters.json"

class arc_length_test(TF.TestFactory):
     file_name = "strategy_tests/arc_length_test/arc_length_test"
     file_parameters = "strategy_tests/arc_length_test/ProjectParameters.json"

class consolidation_interface_3D(TF.TestFactory):
     file_name = "element_tests/consolidation_interface_3D/consolidation_interface_3D"
     file_parameters = "element_tests/consolidation_interface_3D/ProjectParameters.json"


def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            undrained_soil_column_2D,
            arc_length_test,
            consolidation_interface_3D
        ])
    )

    return night_suite
