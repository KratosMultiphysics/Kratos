# Definition of the classes for the VALIDATION TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class rigid_tool_cutting_2D(TF.TestFactory):
    file_name = "solid_tests/cutting_2D/rigid_tool_cutting_2D"
    file_parameters = None

class deformable_tool_cutting_2D(TF.TestFactory):
    file_name = "solid_tests/cutting_2D/deformable_tool_cutting_2D"
    file_parameters = None

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            rigid_tool_cutting_2D,
            deformable_tool_cutting_2D
        ])
    )

    return small_suite
