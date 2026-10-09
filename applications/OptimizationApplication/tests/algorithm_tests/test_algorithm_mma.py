import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_mma import AlgorithmMMA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl


def _CreateModelAndComponents():
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

    properties = model_part.CreateNewProperties(1)
    properties[Kratos.DENSITY] = 2.0
    properties[Kratos.THICKNESS] = 3.0
    model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

    properties = model_part.CreateNewProperties(2)
    properties[Kratos.DENSITY] = 4.0
    properties[Kratos.THICKNESS] = 6.0
    model_part.CreateNewElement("Element2D3N", 2, [4, 1, 3], properties)

    response_function = MassResponseFunction("mass", model, Kratos.Parameters("""{
        "evaluated_model_part_names": ["test"]
    }"""))

    properties_control = MaterialPropertiesControl("control1", model, Kratos.Parameters("""{
        "model_part_names"      : ["test"],
        "control_variable_name" : "DENSITY"
    }"""))
    properties_control.Initialize()

    optimization_problem = OptimizationProblem()
    optimization_problem.AddComponent(response_function)
    optimization_problem.AddComponent(properties_control)
    optimization_problem.AddProcessType("output_processes")

    return model, optimization_problem


class TestAlgorithmMMAUnconstrained(kratos_unittest.TestCase):
    """Minimizes mass w.r.t. DENSITY with no constraints: the box bound alone
    should drive both densities down to their lower bound. Exercises the
    m=0 (unconstrained) short-circuit path in mma_math.solve_mma_subproblem
    through the full Factory/JSON wiring."""

    @classmethod
    def _RunVariant(cls, variant):
        model, optimization_problem = _CreateModelAndComponents()
        parameters = Kratos.Parameters("""{
            "module"    : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"      : "algorithm_mma",
            "objective" : {
                "response_name": "mass",
                "type"         : "minimization",
                "scaling"      : 1.0
            },
            "controls"  : ["control1"],
            "echo_level": 0,
            "settings"  : {
                "echo_level"          : 0,
                "variant"             : \"""" + variant + """\",
                "controls_lower_bound": "0.1",
                "controls_upper_bound": "10.0",
                "conv_settings"       : {
                    "max_iter": 5
                }
            }
        }""")
        algorithm = AlgorithmMMA(model, parameters, optimization_problem)
        algorithm.Initialize()
        conv = algorithm.Solve()
        algorithm.Finalize()
        return algorithm, conv

    def test_mma(self):
        algorithm, conv = self._RunVariant("mma")
        self.assertTrue(conv)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 0.45, places=6)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1, 0.1], places=6)

    def test_gcmma(self):
        algorithm, conv = self._RunVariant("gcmma")
        self.assertTrue(conv)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 0.45, places=6)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1, 0.1], places=6)


class TestAlgorithmMMAConstrained(kratos_unittest.TestCase):
    """Minimizes mass subject to mass >= 5.0 (using the same response as both
    objective and constraint, with a fixed reference value below the
    unconstrained optimum so the constraint is immediately active).
    Exercises the m=1 primal-dual dual-Newton subproblem path through the
    full Factory/JSON wiring, for both the "mma" and "gcmma" variants."""

    @classmethod
    def _RunVariant(cls, variant):
        model, optimization_problem = _CreateModelAndComponents()
        parameters = Kratos.Parameters("""{
            "module"    : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"      : "algorithm_mma",
            "objective" : {
                "response_name": "mass",
                "type"         : "minimization",
                "scaling"      : 1.0
            },
            "constraints": [
                {
                    "response_name"   : "mass",
                    "type"            : ">=",
                    "scaling"         : 1.0,
                    "scaled_ref_value": 5.0
                }
            ],
            "controls"  : ["control1"],
            "echo_level": 0,
            "settings"  : {
                "echo_level"          : 0,
                "variant"             : \"""" + variant + """\",
                "controls_lower_bound": "0.1",
                "controls_upper_bound": "10.0",
                "conv_settings"       : {
                    "max_iter": 10
                }
            }
        }""")
        algorithm = AlgorithmMMA(model, parameters, optimization_problem)
        algorithm.Initialize()
        conv = algorithm.Solve()
        algorithm.Finalize()
        return algorithm, conv

    def test_mma(self):
        algorithm, conv = self._RunVariant("mma")
        self.assertTrue(conv)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 5.0, places=4)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1015, 1.6159], places=3)

    def test_gcmma(self):
        algorithm, conv = self._RunVariant("gcmma")
        self.assertTrue(conv)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 5.0, places=4)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1013, 1.6160], places=3)


if __name__ == "__main__":
    kratos_unittest.main()
