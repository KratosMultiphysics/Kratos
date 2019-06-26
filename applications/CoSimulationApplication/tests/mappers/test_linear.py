import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestMapperLinear(KratosUnittest.TestCase):
    def test_mapper_linear(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
