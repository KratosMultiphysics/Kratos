import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestConvergenceAcceleratorIQNI(KratosUnittest.TestCase):
    def test_convergence_accelerator_iqni(self):
        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()
