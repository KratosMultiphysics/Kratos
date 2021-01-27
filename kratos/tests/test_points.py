import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestPoints(KratosUnittest.TestCase):

    def test_point_constructor_with_kratos_array(self):
        parameters = KratosMultiphysics.Parameters("""{
            "coordinates" : [1.0, 2.0, 3.0]
        }""")
        point = KratosMultiphysics.Point(parameters["coordinates"].GetVector())
        self.assertEqual(point.X, 1.0)
        self.assertEqual(point.Y, 2.0)
        self.assertEqual(point.Z, 3.0)

if __name__ == '__main__':
    KratosUnittest.main()
