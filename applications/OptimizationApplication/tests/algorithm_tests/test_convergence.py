import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import *
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class TestConvergence(kratos_unittest.TestCase):
    @classmethod
    def setUp(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.CreateElements()

        cls.optimization_problem = OptimizationProblem()
        ComponentDataView("algorithm", cls.optimization_problem).SetDataBuffer(10)

    @classmethod
    def CreateElements(cls):
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        properties = cls.model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        properties[Kratos.THICKNESS] = 3.0
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        properties = cls.model_part.CreateNewProperties(2)
        properties[Kratos.DENSITY] = 4.0
        properties[Kratos.THICKNESS] = 6.0
        cls.model_part.CreateNewElement("Element2D3N", 2, [4, 1, 3], properties)

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

    def test_L2(self):
        param = Kratos.Parameters("""{
            "type"              : "l2_norm",
            "max_iter"          : 2,
            "tolerance"         : 1e-3
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        search_direction = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        KratosOA.CollectiveExpressionIO.Read(search_direction, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction / 100000
        self.assertTrue(convergence_criterium.IsConverged())

    def test_L2_max(self):
        param = Kratos.Parameters("""{
            "type"              : "l2_norm",
            "max_iter"          : 2,
            "tolerance"         : 1e-4
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        search_direction = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        KratosOA.CollectiveExpressionIO.Read(search_direction, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertTrue(convergence_criterium.IsConverged())

    def test_AverAbsDelta(self):
        param = Kratos.Parameters("""{
            "type"              : "aver_abs_delta",
            "max_iter"          : 10,
            "tolerance"         : 2e-3
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        algorithm_data.GetBufferedData()["std_obj_value"] = 5
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2222
        self.assertTrue(convergence_criterium.IsConverged())

    def test_AverAbsDelta_max(self):
        param = Kratos.Parameters("""{
            "type"              : "aver_abs_delta",
            "max_iter"          : 2,
            "tolerance"         : 2e-3
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        algorithm_data.GetBufferedData()["std_obj_value"] = 5
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertTrue(convergence_criterium.IsConverged())

    def test_TargetValue(self):
        param = Kratos.Parameters("""{
            "type"              : "target_value",
            "max_iter"          : 10,
            "target_value"      : 0.001
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        algorithm_data.GetBufferedData()["std_obj_value"] = 1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.01
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.001
        self.assertTrue(convergence_criterium.IsConverged())

    def test_TargetValue_max(self):
        param = Kratos.Parameters("""{
            "type"              : "target_value",
            "max_iter"          : 2,
            "target_value"      : 0.001
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = CreateConvergenceCriteria(param, self.optimization_problem)
        algorithm_data.GetBufferedData()["std_obj_value"] = 1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.01
        self.assertTrue(convergence_criterium.IsConverged())

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()