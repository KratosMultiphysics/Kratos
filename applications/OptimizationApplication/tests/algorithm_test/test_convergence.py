from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import *

class TestConvergence(kratos_unittest.TestCase, ABC):
    def test_MaxIter(self):
        param = Kratos.Parameters("""{
            "type"              : "max_iter",
            "max_iter"          : 3
        }""")
        convergence_criterium = CreateConvergenceCriteria(param)
        self.assertFalse(convergence_criterium.CheckConvergence(2))
        self.assertTrue(convergence_criterium.CheckConvergence(3))
        self.assertTrue(convergence_criterium.CheckConvergence(4))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()