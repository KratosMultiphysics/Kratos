import os

# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


# This utility will control the execution scope
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class CoreTests(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)

    def test_array_parameter(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self._test_array_parameter()

    def _test_array_parameter(self):

        default_parameters = KratosMultiphysics.Parameters("""
        {
           "array": [{"one":1, "parameter":{ "one": 1, "array":[1,2], "parameter": { "two": 2 } } },{"two":2},{"three":3}]
        }
        """)
        custom_parameters = KratosMultiphysics.Parameters("""{ }""")

        custom_parameters.AddEmptyValue("array").SetVector(KratosMultiphysics.Vector())

        size = default_parameters["array"].size()
        for i in range(0,size):
            item = default_parameters["array"][i]
            custom_parameters["array"].PushBack(item)

        custom_parameters.ValidateAndAssignDefaults(default_parameters)
        print(" custom_parameters ", custom_parameters.PrettyPrintJsonString())

        item = custom_parameters["array"][size-1]
        if( item.Has("three") == False ):
            raise Exception(" Fail ")


def SetTestSuite(suites):
    core_suite = suites['all']

    core_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CoreTests]))

    return core_suite
