import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine

class TestResponseRoutine(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part_1 = cls.model.CreateModelPart("test1")
        cls.model_part_2 = cls.model.CreateModelPart("test2")
        cls.model_part_3 = cls.model.CreateModelPart("test3")

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test1"],
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control_1 = MaterialPropertiesControl("control1", cls.model, parameters)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test2"],
            "control_variable_name" : "YOUNG_MODULUS"
        }""")
        cls.properties_control_2 = MaterialPropertiesControl("control2", cls.model, parameters)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test1"],
            "control_variable_name" : "THICKNESS"
        }""")
        cls.properties_control_3 = MaterialPropertiesControl("control3", cls.model, parameters)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test3"],
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control_4 = MaterialPropertiesControl("control4", cls.model, parameters)

        for model_part in [cls.model_part_1, cls.model_part_2, cls.model_part_3]:
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            model_part.CreateNewNode(3, 4.0, 4.0, 0.0)

            for i in range(2):
                node_ids = [(i % 3) + 1, ((i + 1) % 3) + 1]
                properties = model_part.CreateNewProperties(i)
                properties[Kratos.DENSITY] = 2.0 * (i + 1)
                properties[Kratos.THICKNESS] = 3.0 * (i + 1)
                properties[Kratos.YOUNG_MODULUS] = 4.0 * (i + 1)
                model_part.CreateNewElement("Element2D2N", i, node_ids, properties)

        cls.master_control = MasterControl()
        cls.master_control.AddControl(cls.properties_control_1)
        cls.master_control.AddControl(cls.properties_control_2)
        cls.master_control.AddControl(cls.properties_control_3)
        cls.master_control.AddControl(cls.properties_control_4)

        cls.properties_control_1.Initialize()
        cls.properties_control_2.Initialize()
        cls.properties_control_3.Initialize()
        cls.properties_control_4.Initialize()

        # now create the response
        parameters = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test1", "test2"
            ]
        }""")
        cls.response = MassResponseFunction("test", cls.model, parameters)
        cls.response.Initialize()

        cls.response_routine = ResponseRoutine(cls.master_control, cls.response)
        cls.response.Initialize()
        cls.response_routine.Initialize()

    def test_CalculateValue(self):
        control_field = self.master_control.GetEmptyField()
        KratosOA.CollectiveExpressionIO.Read(control_field, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change the control field where response does not depend on
        # changing a variable such as YOUNG_MODULUS
        Kratos.Expression.LiteralExpressionIO.SetData(control_field.GetContainerExpressions()[1], 2.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change a dependent variable where the domain is not having intersection
        # changing DENSITY variable
        Kratos.Expression.LiteralExpressionIO.SetData(control_field.GetContainerExpressions()[3], 3.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change a dependent field
        Kratos.Expression.LiteralExpressionIO.SetData(control_field.GetContainerExpressions()[0], 3.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 66)

    def test_CalculateGradient(self):
        control_field = self.master_control.GetEmptyField()
        KratosOA.CollectiveExpressionIO.Read(control_field, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))

        # Calculate value should always be called once before the calculate gradient
        _ = self.response_routine.CalculateValue(control_field)
        gradient = self.response_routine.CalculateGradient()

        self.assertEqual(len(gradient.GetContainerExpressions()), 4)

        # mass response has gradients w.r.t. control1.
        control_1_gradient = gradient.GetContainerExpressions()[0]
        self.assertEqual(Kratos.Expression.Utils.NormL2(control_1_gradient), 20.09975124224178)

        # mass response has gradients w.r.t. YOUNG_MODULUS even the response evaluated in the same control2 domain (same model part).
        control_2_gradient = gradient.GetContainerExpressions()[1]
        self.assertEqual(Kratos.Expression.Utils.NormInf(control_2_gradient), 0.0)

        # mass response has gradients w.r.t. control3.
        control_3_gradient = gradient.GetContainerExpressions()[2]
        self.assertEqual(Kratos.Expression.Utils.NormL2(control_3_gradient), 20.09975124224178)

        # mass response has gradients w.r.t. DENSITY because evaluation and control domains does not have an intersection.
        control_4_gradient = gradient.GetContainerExpressions()[3]
        self.assertEqual(Kratos.Expression.Utils.NormInf(control_4_gradient), 0.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()