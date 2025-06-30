import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.combined_response_function import CombinedResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.geometric_centroid_deviation_response_function import GeometricCentroidDeviationResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction

class TestCombinedResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("../model_part_utils_test/quads_surface", cls.model_part)

        for element in cls.model_part.Elements:
            properties = cls.model_part.CreateNewProperties(element.Id + 5)
            properties[Kratos.DENSITY] = 3.4
            element.Properties = properties

        cls.optimization_problem = OptimizationProblem()

        resp_1_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test.sensitivity_element_1"
            ]
        }""")
        cls.resp_1 = MassResponseFunction("resp_1", cls.model, resp_1_params)
        resp_2_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test.sensitivity_element_2"
            ]
        }""")
        cls.resp_2 = MassResponseFunction("resp_2", cls.model, resp_2_params)
        resp_3_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test.sensitivity_element_3"
            ]
        }""")
        cls.resp_3 = MassResponseFunction("resp_3", cls.model, resp_3_params)
        resp_4_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test.sensitivity_element_4"
            ]
        }""")
        cls.resp_4 = MassResponseFunction("resp_4", cls.model, resp_4_params)
        resp_5_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test.sensitivity_element_2"
            ]
        }""")
        cls.resp_5 = GeometricCentroidDeviationResponseFunction("resp_5", cls.model, resp_5_params)

        cls.optimization_problem.AddComponent(cls.resp_1)
        cls.optimization_problem.AddComponent(cls.resp_2)
        cls.optimization_problem.AddComponent(cls.resp_3)
        cls.optimization_problem.AddComponent(cls.resp_4)
        cls.optimization_problem.AddComponent(cls.resp_5)

        resp_6_params = Kratos.Parameters("""{
            "combining_method"   : "sum",
            "combining_responses": [
                {
                    "response_name"  : "resp_1",
                    "response_weight": 2.1
                },
                {
                    "response_name"  : "resp_2",
                    "response_weight": 3.1
                },
                {
                    "response_name"  : "resp_3",
                    "response_weight": 4.1
                }
            ]
        }""")
        cls.resp_6 = CombinedResponseFunction("resp_6", cls.model, resp_6_params, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.resp_6)

        resp_7_params = Kratos.Parameters("""{
            "combining_method"   : "sum",
            "combining_responses": [
                {
                    "response_name"  : "resp_4",
                    "response_weight": 5.1
                },
                {
                    "response_name"  : "resp_5",
                    "response_weight": 6.1
                }
            ]
        }""")
        cls.resp_7 = CombinedResponseFunction("resp_7", cls.model, resp_7_params, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.resp_7)

        resp_8_params = Kratos.Parameters("""{
            "combining_method"   : "sum",
            "combining_responses": [
                {
                    "response_name"  : "resp_6",
                    "response_weight": 7.1
                },
                {
                    "response_name"  : "resp_7",
                    "response_weight": 8.1
                }
            ]
        }""")
        cls.resp_8 = CombinedResponseFunction("resp_8", cls.model, resp_8_params, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.resp_8)

        cls.resp_8.Initialize()

        cls.model_part.GetNode(1).X = -10.0
        cls.model_part.GetNode(1).Y = -20.0

    def test_CalculateValue(self):
        v1 = self.resp_1.CalculateValue() * 2.1
        v2 = self.resp_2.CalculateValue() * 3.1
        v3 = self.resp_3.CalculateValue() * 4.1
        v4 = self.resp_4.CalculateValue() * 5.1
        v5 = self.resp_5.CalculateValue() * 6.1

        v6 = self.resp_6.CalculateValue()
        v7 = self.resp_7.CalculateValue()

        self.assertAlmostEqual(v6, v1 + v2 + v3)
        self.assertAlmostEqual(v7, v4 + v5)
        self.assertAlmostEqual(self.resp_8.CalculateValue(), v6 * 7.1 + v7 * 8.1)

    def test_CalculateGradient(self):
        def get_data(add_density_physical_var = True):
            elem_exp = Kratos.Expression.ElementExpression(self.model_part)
            nodal_exp = Kratos.Expression.NodalExpression(self.model_part)
            if add_density_physical_var:
                data = {
                    Kratos.DENSITY: KratosOA.CollectiveExpression([elem_exp]),
                    KratosOA.SHAPE: KratosOA.CollectiveExpression([nodal_exp])
                }
            else:
                data = {
                    KratosOA.SHAPE: KratosOA.CollectiveExpression([nodal_exp])
                }
            return nodal_exp, elem_exp, data

        # calculate all the values first
        self.resp_8.CalculateValue()

        v1_shape, v1_rho, data = get_data()
        self.resp_1.CalculateGradient(data)
        v2_shape, v2_rho, data = get_data()
        self.resp_2.CalculateGradient(data)
        v3_shape, v3_rho, data = get_data()
        self.resp_3.CalculateGradient(data)
        v4_shape, v4_rho, data = get_data()
        self.resp_4.CalculateGradient(data)
        v5_shape, v5_rho, data = get_data(False)
        self.resp_5.CalculateGradient(data)

        v6_shape, v6_rho, data = get_data()
        self.resp_6.CalculateGradient(data)
        v7_shape, v7_rho, data = get_data()
        self.resp_7.CalculateGradient(data)

        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v6_shape - v1_shape * 2.1 - v2_shape * 3.1 - v3_shape * 4.1), 0.0)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v6_rho - v1_rho * 2.1 - v2_rho * 3.1 - v3_rho * 4.1), 0.0)

        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v7_shape - v4_shape * 5.1 - v5_shape * 6.1), 0.0)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v7_rho - v4_rho * 5.1), 0.0)

        v8_shape, v8_rho, data = get_data()
        self.resp_8.CalculateGradient(data)

        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v8_shape - v6_shape * 7.1 - v7_shape * 8.1), 0.0)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(v8_rho - v6_rho * 7.1 - v7_rho * 8.1), 0.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()