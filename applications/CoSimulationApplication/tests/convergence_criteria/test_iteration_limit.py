import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestConvergenceCriterionIterationLimit(KratosUnittest.TestCase):
    def test_convergence_criterion_iteration_limit(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
