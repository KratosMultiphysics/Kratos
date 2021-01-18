# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class CoreTests(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)

    @classmethod
    def test_array_parameter(self):

        default_parameters = KratosMultiphysics.Parameters("""
        {
           "array": [{"one":1, "parameter":{ "one": 1, "array":[1,2], "parameter": { "two": 2 } } },{"two":2},{"three":3}]
        }
        """)
        custom_parameters = KratosMultiphysics.Parameters("""{ }""")

        custom_parameters.AddEmptyList("array")

        size = default_parameters["array"].size()
        for i in range(0,size):
            item = default_parameters["array"][i]
            custom_parameters["array"].Append(item)

        #append double
        custom_parameters["array"].Append(500.1)
        #append int
        custom_parameters["array"].Append(80)
        #append int
        custom_parameters["array"].Append(0)
        #append bool
        custom_parameters["array"].Append(True)
        #append string
        custom_parameters["array"].Append("hello")
        #append vector
        custom_parameters["array"].Append(KratosMultiphysics.Vector(2))
        #append matrix
        custom_parameters["array"].Append(KratosMultiphysics.Matrix(2,2))

        #custom_parameters.ValidateAndAssignDefaults(default_parameters)
        print(" custom_parameters ", custom_parameters.PrettyPrintJsonString())

        item = custom_parameters["array"][size-1]
        if( item.Has("three") == False ):
            raise Exception(" Fail ")


def SetTestSuite(suites):
    core_suite = suites['all']

    core_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CoreTests]))

    return core_suite
