import numpy as np

import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl

# The GP algorithms import LinearSolversApplication at module level, so they are imported only when it is available.
if CheckIfApplicationsAvailable("LinearSolversApplication"):
    from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_gradient_projection import AlgorithmGradientProjection
    from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_relaxed_gradient_projection import AlgorithmRelaxedGradientProjection


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
    """Objective and constraint gradients are parallel (both a mass), start at mass = 15.
    Violated mass <= 1: GP's projected direction vanishes, so kkt_norm has to equal the violation; RGP adds the
    correction term to the direction. Active mass <= 15 / mass >= 15: kkt_norm must see the sign of the inequality
    multiplier (dual feasibility). GP: kkt_norm = max(|g - N lambda_clipped|_inf, violation); RGP:
    kkt_norm = max(|search_direction|_inf, violation)."""

    def _Run(self, algorithm_class, extra_settings="", constraint_settings="", constraint_type="<=", reference=1.0, check_violated=True):
        model, optimization_problem = _CreateModelAndProblem()
        parameters = Kratos.Parameters("""{
            "module"     : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"       : "gradient_projection",
            "objective"  : {"response_expression": "mass", "type": "minimization", "scaling": 1.0},
            "constraints": [{"response_expression": "mass_copy", "type": \"""" + constraint_type + """\", "scaled_ref_value": """ + str(reference) + constraint_settings + """}],
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
        if check_violated:
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

    def _KKTNormAtActiveConstraint(self, algorithm_class, constraint_type):
        # mass = 15 at the start, so mass_copy <= / >= 15 is active (not violated); grad(mass) = [1.5, 3], parallel to the constraint gradient.
        kkt_norm, _, violation = self._Run(algorithm_class, constraint_type=constraint_type, reference=15.0, check_violated=False)
        self.assertLess(abs(violation), 1e-12)
        return kkt_norm

    def test_gradient_projection_dual_infeasible(self):
        # Minimising mass with mass <= 15 active: the multiplier has the wrong sign (the constraint should be released),
        # the projected direction is 0, but the point is not KKT: kkt_norm = |grad|_inf = 3.
        self.assertAlmostEqual(self._KKTNormAtActiveConstraint(AlgorithmGradientProjection, "<="), 3.0, 12)

    def test_gradient_projection_kkt_point(self):
        # With mass >= 15 the multiplier is admissible and the stationarity residual vanishes.
        self.assertLess(self._KKTNormAtActiveConstraint(AlgorithmGradientProjection, ">="), 1e-12)

    def test_relaxed_gradient_projection_dual_infeasible(self):
        # RGP already clips the multiplier in its projection (the direction is non-zero), so kkt_norm is not 0 either.
        self.assertAlmostEqual(self._KKTNormAtActiveConstraint(AlgorithmRelaxedGradientProjection, "<="), 1.0, 12)


if __name__ == "__main__":
    kratos_unittest.main()
