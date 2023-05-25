import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import *
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class TestConvergence(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.optimization_problem = OptimizationProblem()

    def test_MaxIter(self):
        param = Kratos.Parameters("""{
            "type"              : "max_iter",
            "max_iter"          : 2
        }""")
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertTrue(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertTrue(convergence_criterium.IsConverged())

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()