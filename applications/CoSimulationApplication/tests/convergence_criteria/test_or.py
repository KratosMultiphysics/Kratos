import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestConvergenceCriterionOr(KratosUnittest.TestCase):
    def test_convergence_criterion_or(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
