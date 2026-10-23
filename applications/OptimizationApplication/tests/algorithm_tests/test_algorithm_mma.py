import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_mma import AlgorithmMMA
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl


def _CreateModelAndComponents():
    """Two triangles in sub model parts part_a / part_b (masses 1.5 rho_1 and 3.0 rho_2, responses mass, mass_A, mass_B), one DENSITY control each."""
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

    model_part.CreateSubModelPart("part_a").AddElements([1])
    model_part.CreateSubModelPart("part_b").AddElements([2])

    optimization_problem = OptimizationProblem()
    for name, part in (("mass", "test"), ("mass_A", "test.part_a"), ("mass_B", "test.part_b")):
        optimization_problem.AddComponent(MassResponseFunction(name, model, Kratos.Parameters("""{
            "evaluated_model_part_names": [\"""" + part + """\"]
        }""")))

    for name, part in (("control_a", "test.part_a"), ("control_b", "test.part_b")):
        control = MaterialPropertiesControl(name, model, Kratos.Parameters("""{
            "model_part_names"      : [\"""" + part + """\"],
            "control_variable_name" : "DENSITY"
        }"""))
        control.Initialize()
        optimization_problem.AddComponent(control)

    optimization_problem.AddProcessType("output_processes")
    return model, optimization_problem


def _Parameters(variant, objective, constraints="", bounds=None, max_iter=20, extra_settings=""):
    if bounds is None:
        bounds = [("control_a", 0.1, 10.0), ("control_b", 0.1, 10.0)]
    bounds_json = ",".join(f'{{"control": "{name}", "lower": {lower}, "upper": {upper}}}' for name, lower, upper in bounds)
    return Kratos.Parameters("""{
        "module"     : "KratosMultiphysics.OptimizationApplication.algorithms",
        "type"       : "algorithm_mma",
        "objective"  : {"response_expression": \"""" + objective + """\", "type": "minimization", "scaling": 1.0},
        "constraints": [""" + constraints + """],
        "controls"   : ["control_a", "control_b"],
        "settings"   : {
            "variant"       : \"""" + variant + """\",
            "control_bounds": [""" + bounds_json + """],
            "conv_settings" : {"max_iter": """ + str(max_iter) + """}""" + extra_settings + """
        }
    }""")


class _FieldRecorder(Kratos.OutputProcess):
    """Output process recording the physical DENSITY of the model next to the algorithm's control field."""
    def __init__(self, model, optimization_problem):
        super().__init__()
        self.model_part = model["test"]
        self.optimization_problem = optimization_problem
        self.records = []

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        data = ComponentDataView("algorithm", self.optimization_problem).GetBufferedData()
        physical = [element.Properties[Kratos.DENSITY] for element in self.model_part.Elements]
        self.records.append((physical, list(data["control_field"].data), data["gcmma_inner_iterations"]))


def _Run(parameters, recorder=False):
    model, optimization_problem = _CreateModelAndComponents()
    if recorder:
        optimization_problem.AddProcess("output_processes", _FieldRecorder(model, optimization_problem))
    algorithm = AlgorithmMMA(model, parameters, optimization_problem)
    algorithm.Initialize()
    converged = algorithm.Solve()
    algorithm.Finalize()
    if recorder:
        return optimization_problem.GetListOfProcesses("output_processes")[0].records
    return algorithm, optimization_problem, converged


class TestAlgorithmMMAUnconstrained(kratos_unittest.TestCase):
    """Minimizing the total mass: each density goes to the lower bound of its own control."""

    def _CheckVariant(self, variant):
        bounds = [("control_a", 0.1, 10.0), ("control_b", 0.5, 8.0)]
        algorithm, optimization_problem, converged = _Run(_Parameters(variant, "mass", bounds=bounds, max_iter=8))
        self.assertTrue(converged)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1, 0.5], 10)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 1.5 * 0.1 + 3.0 * 0.5, 10)
        # both variables are at their bounds, so the KKT residual has to vanish
        data = ComponentDataView("algorithm", optimization_problem).GetBufferedData()
        self.assertLess(data["kkt_norm"], 1e-9)
        self.assertEqual(data["artificial_variable_max"], 0.0)

    def test_mma(self):
        self._CheckVariant("mma")

    def test_gcmma(self):
        self._CheckVariant("gcmma")

    def test_upper_bounds(self):
        # maximising 1/mass is not available, so drive both densities up through a violated constraint
        constraint = '{"response_expression": "mass", "type": ">=", "scaling": 1.0, "scaled_ref_value": 1000.0}'
        bounds = [("control_a", 0.1, 3.0), ("control_b", 0.1, 4.0)]
        algorithm, _, _ = _Run(_Parameters("mma", "mass", constraint, bounds=bounds, max_iter=12, extra_settings=', "move": 0.9'))
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [3.0, 4.0], 10)


class TestAlgorithmMMAConstrained(kratos_unittest.TestCase):
    """Constrained problems with a unique optimum (an objective that ignores a control would leave it undetermined)."""

    def _CheckBoundAndActiveConstraint(self, variant):
        # min mass s.t. mass_A >= 0.9: rho_1 = 0.9 / 1.5 (active constraint), rho_2 = 0.1 (lower bound)
        constraint = '{"response_expression": "mass_A", "type": ">=", "scaling": 1.0, "scaled_ref_value": 0.9}'
        algorithm, optimization_problem, converged = _Run(_Parameters(variant, "mass", constraint, max_iter=30))
        self.assertTrue(converged)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.6, 0.1], 6)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 1.5 * 0.6 + 3.0 * 0.1, 6)
        data = ComponentDataView("algorithm", optimization_problem).GetBufferedData()
        self.assertLess(data["kkt_norm"], 1e-5)
        self.assertEqual(data["artificial_variable_max"], 0.0)

    def _CheckInteriorOptimum(self, variant):
        # min mass s.t. mass_A * mass_B >= 9, i.e. rho_1 rho_2 >= 2: optimum rho = (2, 1), both away from their bounds
        constraint = '{"response_expression": "mass_A * mass_B", "type": ">=", "scaling": 1.0, "scaled_ref_value": 9.0}'
        algorithm, optimization_problem, converged = _Run(_Parameters(variant, "mass", constraint, max_iter=40, extra_settings=', "move": 0.3'))
        self.assertTrue(converged)
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [2.0, 1.0], 5)
        self.assertAlmostEqual(algorithm.GetOptimizedObjectiveValue(), 6.0, 5)
        self.assertLess(ComponentDataView("algorithm", optimization_problem).GetBufferedData()["kkt_norm"], 1e-4)

    def test_mma_bound_and_active_constraint(self):
        self._CheckBoundAndActiveConstraint("mma")

    def test_gcmma_bound_and_active_constraint(self):
        self._CheckBoundAndActiveConstraint("gcmma")

    def test_mma_interior_optimum(self):
        self._CheckInteriorOptimum("mma")

    def test_gcmma_interior_optimum(self):
        self._CheckInteriorOptimum("gcmma")

    def test_gcmma_output_fields_match_control_field(self):
        # Trial evaluations leave the physical fields at the trial design; the output of step k must show x_k.
        constraint = '{"response_expression": "mass_A * mass_B", "type": ">=", "scaling": 1.0, "scaled_ref_value": 9.0}'
        records = _Run(_Parameters("gcmma", "mass", constraint, max_iter=6, extra_settings=', "move": 0.9'), recorder=True)
        self.assertGreater(len(records), 3)
        self.assertTrue(any(inner > 0 for _, _, inner in records))   # otherwise the test cannot detect a stale field
        for physical, control_field, _ in records:
            self.assertVectorAlmostEqual(physical, control_field, 12)

    def test_zero_gradient_control_does_not_move(self):
        # mass_A does not depend on control_b: with the regularisation of p, q it must stay at its initial value.
        constraint = '{"response_expression": "mass", "type": ">=", "scaling": 1.0, "scaled_ref_value": 5.0}'
        algorithm, _, _ = _Run(_Parameters("mma", "mass_A", constraint, max_iter=10))
        self.assertVectorAlmostEqual(algorithm.GetCurrentControlField().data, [0.1, 4.0], 10)

    def test_infeasible_constraint_activates_artificial_variable(self):
        constraint = '{"response_expression": "mass", "type": ">=", "scaling": 1.0, "scaled_ref_value": 1e6}'
        _, optimization_problem, _ = _Run(_Parameters("mma", "mass_A", constraint, max_iter=2))
        data = ComponentDataView("algorithm", optimization_problem).GetBufferedData()
        self.assertGreater(data["artificial_variable_max"], 1.0)
        self.assertGreater(data["kkt_norm"], 1.0)


class TestAlgorithmMMAValidation(kratos_unittest.TestCase):
    def _Create(self, parameters):
        model, optimization_problem = _CreateModelAndComponents()
        return AlgorithmMMA(model, parameters, optimization_problem), model, optimization_problem

    def test_missing_control_bounds_entry(self):
        with self.assertRaisesRegex(RuntimeError, "requires a \"control_bounds\" entry for control \"control_b\""):
            self._Create(_Parameters("mma", "mass", bounds=[("control_a", 0.1, 1.0)]))

    def test_unknown_control_in_bounds(self):
        with self.assertRaisesRegex(RuntimeError, "unknown control \"other\""):
            self._Create(_Parameters("mma", "mass", bounds=[("control_a", 0.1, 1.0), ("control_b", 0.1, 1.0), ("other", 0.1, 1.0)]))

    def test_lower_not_below_upper(self):
        with self.assertRaisesRegex(RuntimeError, "requires lower < upper"):
            self._Create(_Parameters("mma", "mass", bounds=[("control_a", 1.0, 1.0), ("control_b", 0.1, 1.0)]))

    def test_incomplete_entry(self):
        parameters = _Parameters("mma", "mass")
        parameters["settings"]["control_bounds"][0].RemoveValue("upper")
        with self.assertRaisesRegex(RuntimeError, "require"):
            self._Create(parameters)

    def test_initial_design_outside_bounds(self):
        algorithm, _, _ = self._Create(_Parameters("mma", "mass", bounds=[("control_a", 0.1, 1.0), ("control_b", 0.1, 10.0)]))
        with self.assertRaisesRegex(RuntimeError, "outside \"control_bounds\""):
            algorithm.Initialize()   # DENSITY of part_a is 2.0 > 1.0

    def test_invalid_variant_and_equality(self):
        with self.assertRaisesRegex(RuntimeError, "Unsupported MMA"):
            self._Create(_Parameters("foo", "mass"))
        equality = '{"response_expression": "mass", "type": "=", "scaling": 1.0, "scaled_ref_value": 5.0}'
        with self.assertRaisesRegex(RuntimeError, "inequality constraints only"):
            self._Create(_Parameters("mma", "mass", equality))

    def test_invalid_artificial_coefficients(self):
        with self.assertRaisesRegex(RuntimeError, "\"d\" > 0"):
            self._Create(_Parameters("mma", "mass", extra_settings=', "d": 0.0'))

    def test_invalid_raa0(self):
        with self.assertRaisesRegex(RuntimeError, "\"raa0\" > 0"):
            self._Create(_Parameters("mma", "mass", extra_settings=', "raa0": 0.0'))


if __name__ == "__main__":
    kratos_unittest.main()
