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
            custom_parameters["array"].__setitem__(i,default_parameters["array"].__getitem__(i))

        if( custom_parameters["array"][2].Has("three") ):
            return True
        else:
            return False



def SetTestSuite(suites):
    core_suite = suites['small']

    core_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CoreTests]))

    return core_suite
