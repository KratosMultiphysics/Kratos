# import Kratos
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class CoreTests(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)


    def test_array_parameter(self):

        default_parameters = KratosMultiphysics.Parameters("""
        {
           "array": [{"one":1},{"two":2},{"three":3}]
        }
        """)
        custom_parameters = KratosMultiphysics.Parameters("""{  }""")
        size = 3
        custom_parameters.AddEmptyValue("array").SetVector(KratosMultiphysics.Vector(size))

        for i in range(0,size):
            item = default_parameters["array"][i]
            custom_parameters["array"][i] = item

        item = custom_parameters["array"][size-1]
        if( item.Has("three") == False ):
            raise Exception(" Fail 2 ")


def SetTestSuite(suites):
    core_suite = suites['small']

    core_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CoreTests]))

    return core_suite
