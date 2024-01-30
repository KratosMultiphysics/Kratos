import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import *
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestLineSearch(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.CreateElements()

        cls.optimization_problem = OptimizationProblem()
        cls.optimization_problem_one_step = OptimizationProblem()
        ComponentDataView("algorithm", cls.optimization_problem).SetDataBuffer(2)
        ComponentDataView("algorithm", cls.optimization_problem_one_step).SetDataBuffer(2)

        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(cls.model_part)])
        KratosOA.CollectiveExpressionIO.Read(sensitivity, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))
        ComponentDataView("algorithm", cls.optimization_problem).GetBufferedData()["search_direction"] = sensitivity * 10
        ComponentDataView("algorithm", cls.optimization_problem_one_step).GetBufferedData()["search_direction"] = sensitivity
        ComponentDataView("algorithm", cls.optimization_problem).GetBufferedData()["control_field_update"] = sensitivity * 0.2
        cls.optimization_problem.AdvanceStep()
        ComponentDataView("algorithm", cls.optimization_problem).GetBufferedData()["search_direction"] = sensitivity
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

    def test_ConstantLineSearchInfNorm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "const_step",
            "gradient_scaling": "inf_norm",
            "init_step"          : 3.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.75)

    def test_ConstantLineSearchL2Norm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "const_step",
            "gradient_scaling": "l2_norm",
            "init_step"          : 3.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.6708203932499369)

    def test_ConstantLineSearchNoneNorm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "const_step",
            "gradient_scaling"  : "none",
            "init_step"         : 3.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 3.0)

    def test_BBStepDefParam(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "inf_norm",
            "init_step"         : 0.0,
            "max_step"          : 0.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.0)

    def test_BBStepInfNorm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "inf_norm",
            "init_step"         : 3.0,
            "max_step"          : 10.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.005555555555555557)

    def test_BBStepL2Norm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "l2_norm",
            "init_step"         : 3.0,
            "max_step"          : 10.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.004969039949999534)

    def test_BBStepNoneNorm(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "none",
            "init_step"         : 3.0,
            "max_step"          : 10.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.022222222222222227)

    def test_BBStepNoneNormMaxStep(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "none",
            "init_step"         : 0.001,
            "max_step"          : 0.01
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 0.01)

    def test_BBStepNoneNormInitStep(self):
        line_search_settings = Kratos.Parameters("""{
            "type"              : "BB_step",
            "gradient_scaling"  : "none",
            "init_step"         : 3.0,
            "max_step"          : 10.0
        }""")
        line_search = CreateLineSearch(line_search_settings, self.optimization_problem_one_step)
        alpha = line_search.ComputeStep()
        self.assertEqual(alpha, 3.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()