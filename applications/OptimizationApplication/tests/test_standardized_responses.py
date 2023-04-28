
from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.responses.standardized_constraint import StandardizedConstraint

class TestStandardizedComponent(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.optimization_info = OptimizationProblem()

        cls.response_function = MassResponseFunction(cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""), cls.optimization_info)
        cls.optimization_info.AddResponse("mass", cls.response_function)

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

        cls.response_function.Initialize()
        cls.response_function.Check()

    def _CheckSensitivity(self, standardized_component: Union[StandardizedObjective, StandardizedConstraint], delta: float, precision: int):
        self.optimization_info.AdvanceStep()
        container_expression = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        collective_expression = KratosOA.ContainerExpression.CollectiveExpressions([container_expression])
        sensitivities = {KratosOA.DENSITY_SENSITIVITY: collective_expression}
        standardized_component.CalculateStandardizedSensitivity(sensitivities)
        container_expression.Evaluate(KratosOA.YOUNG_MODULUS_SENSITIVITY)

        ref_value = standardized_component.GetStandardizedValue()
        for element in self.model_part.Elements:
            element.Properties[Kratos.DENSITY] += delta
            self.optimization_info.AdvanceStep()
            value = standardized_component.GetStandardizedValue()
            sensitivity = (value - ref_value)/delta
            element.Properties[Kratos.DENSITY] -= delta

            ref_sensitivity = element.Properties[KratosOA.YOUNG_MODULUS_SENSITIVITY]
            self.assertAlmostEqual(ref_sensitivity, sensitivity, precision)

class TestStandardizedObjective(TestStandardizedComponent):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "maximization",
            "scaling"      : 2.0
        }""")

        cls.standardized_objective = StandardizedObjective(parameters, cls.optimization_info)

    def test_ObjectiveMinimization(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "minimization",
            "scaling"      : 2.0
        }""")

        standardized_objective = StandardizedObjective(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_objective.GetStandardizedValue(), self.response_function.CalculateValue() * 2.0)
        self._CheckSensitivity(standardized_objective, 1e-9, 6)

    def test_ObjectiveMaximization(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "maximization",
            "scaling"      : 2.0
        }""")

        standardized_objective = StandardizedObjective(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_objective.GetStandardizedValue(), self.response_function.CalculateValue() * -2.0)
        self._CheckSensitivity(standardized_objective, 1e-9, 6)

    def test_GetInitialValue(self):
        self.assertEqual(self.standardized_objective.GetInitialValue(), self.response_function.CalculateValue())

    def test_GetResponseFunctionName(self):
        self.assertEqual(self.standardized_objective.GetName(), "mass")

    def test_GetResponseType(self):
        self.assertEqual(self.standardized_objective.GetType(), "maximization")

    def test_UpdateObjectiveData(self):
        self.standardized_objective.UpdateObjectiveData()
        response_data: BufferedDict = self.optimization_info.GetReponseData("mass")["buffered"]
        self.assertEqual(response_data["type"], self.standardized_objective.GetType())
        self.assertTrue(response_data.HasValue("rel_change"))
        self.assertTrue(response_data.HasValue("abs_change"))

class TestStandardizedConstraint(TestStandardizedComponent):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : ">=",
            "ref_value"    : 4.0
        }""")

        cls.standardized_constraint = StandardizedConstraint(parameters, cls.optimization_info)

    def test_ConstraintInitialValueEquality(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "=",
            "ref_value"    : "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self.assertTrue(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintInitialValueLessThan(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "<=",
            "ref_value"    : "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self.assertTrue(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintInitialValueGreaterThan(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : ">=",
            "ref_value"    : "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self.assertTrue(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueEquality(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "=",
            "ref_value"    : 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), self.response_function.CalculateValue() - 4.0)
        self.assertTrue(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueLessThan(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "<=",
            "ref_value"    : 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), self.response_function.CalculateValue() - 4.0)
        self.assertTrue(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueGreaterThan(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : ">=",
            "ref_value"    : 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.optimization_info)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), -self.response_function.CalculateValue() + 4.0)
        self.assertFalse(standardized_constraint.IsActive())
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_GetReferenceValue(self):
        self.assertEqual(self.standardized_constraint.GetReferenceValue(), 4.0)

    def test_GetResponseFunctionName(self):
        self.assertEqual(self.standardized_constraint.GetName(), "mass")

    def test_GetResponseType(self):
        self.assertEqual(self.standardized_constraint.GetType(), ">=")

    def test_UpdateConstraintData(self):
        self.standardized_constraint.UpdateConstraintData()
        response_data: BufferedDict = self.optimization_info.GetReponseData("mass")["buffered"]
        self.assertEqual(response_data["type"], self.standardized_constraint.GetType())
        self.assertTrue(response_data.HasValue("rel_change"))
        self.assertTrue(response_data.HasValue("abs_change"))
        self.assertTrue(response_data.HasValue("violation"))
        self.assertFalse(response_data["is_active"])
        self.assertEqual(response_data["ref_value"], 4.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()