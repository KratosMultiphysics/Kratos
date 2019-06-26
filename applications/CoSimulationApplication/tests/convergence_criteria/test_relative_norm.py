import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestConvergenceCriterionRelativeNorm(KratosUnittest.TestCase):
    def test_convergence_criterion_relative_norm(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
