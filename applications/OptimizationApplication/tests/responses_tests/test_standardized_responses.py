
from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestStandardizedComponent(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.optimization_problem = OptimizationProblem()

        cls.response_function = MassResponseFunction("mass", cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""))
        cls.optimization_problem.AddComponent(cls.response_function)

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

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test"],
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control = MaterialPropertiesControl("control1", cls.model, parameters)

        cls.master_control = MasterControl()
        cls.master_control.AddControl(cls.properties_control)
        cls.properties_control.Initialize()

        cls.initial_configuration = cls.master_control.GetEmptyField()
        KratosOA.CollectiveExpressionIO.Read(cls.initial_configuration, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))

    def _CheckSensitivity(self, standardized_component: Union[StandardizedObjective, StandardizedConstraint], delta: float, precision: int):
        self.optimization_problem.AdvanceStep()
        ref_value = standardized_component.CalculateStandardizedValue(self.initial_configuration)
        gradients = standardized_component.CalculateStandardizedGradient()
        KratosOA.CollectiveExpressionIO.Write(gradients, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.YOUNG_MODULUS))

        for element in self.model_part.Elements:
            element.Properties[Kratos.DENSITY] += delta
            current_configuration = self.master_control.GetControlField()
            value = standardized_component.CalculateStandardizedValue(current_configuration, False)
            sensitivity = (value - ref_value)/delta
            element.Properties[Kratos.DENSITY] -= delta

            ref_sensitivity = element.Properties[Kratos.YOUNG_MODULUS]
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

        cls.standardized_objective = StandardizedObjective(parameters, cls.master_control, cls.optimization_problem)
        cls.standardized_objective.CalculateStandardizedValue(cls.initial_configuration)

    def test_ObjectiveMinimization(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "minimization",
            "scaling"      : 2.0
        }""")

        standardized_objective = StandardizedObjective(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_objective.GetStandardizedValue(), self.response_function.CalculateValue() * 2.0)
        self._CheckSensitivity(standardized_objective, 1e-9, 6)

    def test_ObjectiveMaximization(self):
        parameters = Kratos.Parameters("""{
            "response_name": "mass",
            "type"         : "maximization",
            "scaling"      : 2.0
        }""")

        standardized_objective = StandardizedObjective(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_objective.GetStandardizedValue(), self.response_function.CalculateValue() * -2.0)
        self._CheckSensitivity(standardized_objective, 1e-9, 6)

    def test_GetInitialValue(self):
        self.assertEqual(self.standardized_objective.GetInitialValue(), self.response_function.CalculateValue() * -2.0)

class TestStandardizedConstraint(TestStandardizedComponent):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : ">=",
            "scaled_ref_value"    : 4.0
        }""")

        cls.standardized_constraint = StandardizedConstraint(parameters, cls.master_control, cls.optimization_problem)
        cls.standardized_constraint.CalculateStandardizedValue(cls.initial_configuration)

    def test_ConstraintInitialValueEquality(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : "=",
            "scaled_ref_value": "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintInitialValueLessThan(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : "<=",
            "scaled_ref_value": "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintInitialValueGreaterThan(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : ">=",
            "scaled_ref_value": "initial_value"
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), 0.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueEquality(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : "=",
            "scaled_ref_value": 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), self.response_function.CalculateValue() - 4.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueLessThan(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : "<=",
            "scaled_ref_value": 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), self.response_function.CalculateValue() - 4.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_ConstraintSpecifiedValueGreaterThan(self):
        parameters = Kratos.Parameters("""{
            "response_name"   : "mass",
            "type"            : ">=",
            "scaled_ref_value": 4.0
        }""")

        standardized_constraint = StandardizedConstraint(parameters, self.master_control, self.optimization_problem)
        self.assertAlmostEqual(standardized_constraint.GetStandardizedValue(), -self.response_function.CalculateValue() + 4.0)
        self._CheckSensitivity(standardized_constraint, 1e-9, 6)

    def test_GetStandardizedReferenceValue(self):
        self.assertEqual(self.standardized_constraint.GetStandardizedReferenceValue(), -4.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()