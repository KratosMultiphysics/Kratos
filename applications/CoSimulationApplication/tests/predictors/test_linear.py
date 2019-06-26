import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPredictorLinear(KratosUnittest.TestCase):
    def test_predictor_linear(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
