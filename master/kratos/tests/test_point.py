import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPoint(KratosUnittest.TestCase):

    def test_point_constructor_with_kratos_array(self):
        coords = [1.0, -2.5, 3.3]
        arr = KM.Array3(coords)

        point = KM.Point(arr)
        self.assertAlmostEqual(point.X, coords[0])
        self.assertAlmostEqual(point.Y, coords[1])
        self.assertAlmostEqual(point.Z, coords[2])

    def test_point_constructor_with_kratos_vector(self):
        coords = [1.0, -2.5, 3.3]
        vec = KM.Vector(coords)

        point = KM.Point(vec)
        self.assertAlmostEqual(point.X, coords[0])
        self.assertAlmostEqual(point.Y, coords[1])
        self.assertAlmostEqual(point.Z, coords[2])


if __name__ == '__main__':
    KratosUnittest.main()
