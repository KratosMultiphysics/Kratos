import numpy as np

import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_gradient_projection import AlgorithmGradientProjection
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_relaxed_gradient_projection import AlgorithmRelaxedGradientProjection
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl


def _CreateModelAndProblem():
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
    for element_id, nodes, density, thickness in ((1, [1, 2, 3], 2.0, 3.0), (2, [4, 1, 3], 4.0, 6.0)):
        properties = model_part.CreateNewProperties(element_id)
        properties[Kratos.DENSITY] = density
        properties[Kratos.THICKNESS] = thickness
        model_part.CreateNewElement("Element2D3N", element_id, nodes, properties)

    optimization_problem = OptimizationProblem()
    for name in ("mass", "mass_copy"):   # two responses so objective and constraint get their own data buffers
        optimization_problem.AddComponent(MassResponseFunction(name, model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}""")))
    control = MaterialPropertiesControl("control1", model, Kratos.Parameters("""{
        "model_part_names": ["test"], "control_variable_name": "DENSITY"}"""))
    control.Initialize()
    optimization_problem.AddComponent(control)
    optimization_problem.AddProcessType("output_processes")
    return model, optimization_problem


@kratos_unittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestGradientProjectionKKTNorm(kratos_unittest.TestCase):
    """Objective and constraint gradients are parallel (both a mass) and the constraint mass <= 1 is violated at the
    start (mass = 15). GP: the projected direction vanishes, so kkt_norm has to equal the violation. RGP: the
    correction term is added to the direction. In both cases kkt_norm = max(|search_direction|_inf, violation)."""

    def _Run(self, algorithm_class, extra_settings="", constraint_settings=""):
        model, optimization_problem = _CreateModelAndProblem()
        parameters = Kratos.Parameters("""{
            "module"     : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"       : "gradient_projection",
            "objective"  : {"response_expression": "mass", "type": "minimization", "scaling": 1.0},
            "constraints": [{"response_expression": "mass_copy", "type": "<=", "scaled_ref_value": 1.0""" + constraint_settings + """}],
            "controls"   : ["control1"],
            "settings"   : {
                "line_search"  : {"type": "const_step", "init_step": 0.1, "gradient_scaling": "inf_norm"},
                "conv_settings": {"max_iter": 1}""" + extra_settings + """
            }
        }""")
        algorithm = algorithm_class(model, parameters, optimization_problem)
        algorithm.Initialize()
        algorithm.Solve()
        algorithm.Finalize()
        data = ComponentDataView("algorithm", optimization_problem).GetBufferedData()
        violation = algorithm._constraints_list[0].GetStandardizedValue()
        self.assertGreater(violation, 10.0)
        search_direction = data["search_direction"].data
        return data["kkt_norm"], search_direction, violation

    def test_gradient_projection(self):
        kkt_norm, search_direction, violation = self._Run(AlgorithmGradientProjection)
        self.assertLess(np.max(np.abs(search_direction)), 1e-12)   # parallel gradients: nothing left to project
        self.assertAlmostEqual(kkt_norm, violation, 12)

    def test_relaxed_gradient_projection(self):
        # a small correction cap keeps the direction below the violation, so the violation has to set kkt_norm
        kkt_norm, search_direction, violation = self._Run(AlgorithmRelaxedGradientProjection, constraint_settings=', "max_w_c": 0.01')
        self.assertLess(np.max(np.abs(search_direction)), violation)
        self.assertAlmostEqual(kkt_norm, violation, 12)
        # the default cap lets the correction dominate: kkt_norm is then the direction norm
        kkt_norm, search_direction, violation = self._Run(AlgorithmRelaxedGradientProjection)
        self.assertAlmostEqual(kkt_norm, max(np.max(np.abs(search_direction)), violation), 12)
        self.assertGreaterEqual(kkt_norm, violation)


if __name__ == "__main__":
    kratos_unittest.main()
