import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
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
        cls.optimization_problem = OptimizationProblem()

        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>",
            "model_part_names"      : ["test1"],
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control_1 = MaterialPropertiesControl(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddControl("control1", cls.properties_control_1)

        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>",
            "model_part_names"      : ["test2"],
            "control_variable_name" : "YOUNG_MODULUS"
        }""")
        cls.properties_control_2 = MaterialPropertiesControl(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddControl("control2", cls.properties_control_2)

        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>",
            "model_part_names"      : ["test1"],
            "control_variable_name" : "THICKNESS"
        }""")
        cls.properties_control_3 = MaterialPropertiesControl(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddControl("control3", cls.properties_control_3)

        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>",
            "model_part_names"      : ["test3"],
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control_4 = MaterialPropertiesControl(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddControl("control4", cls.properties_control_4)

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

        cls.master_control.Initialize()

        # now create the response
        parameters = Kratos.Parameters("""{
            "combined_output_model_part_name": "<RESPONSE_NAME>",
            "evaluated_model_part_names"     : [
                "test1", "test2"
            ]
        }""")
        cls.response = MassResponseFunction(cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddResponse("test", cls.response)

        cls.response_routine = ResponseRoutine(cls.master_control, "test", cls.optimization_problem)
        cls.response_routine.Initialize()

    def test_CalculateValue(self):
        control_field = self.master_control.GetEmptyControlFields()
        control_field.Read(Kratos.DENSITY)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change the control field where response does not depend on
        # changing a variable such as YOUNG_MODULUS
        control_field.GetContainerExpressions()[1].SetData(2.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change a dependent variable where the domain is not having intersection
        # changing DENSITY variable
        control_field.GetContainerExpressions()[3].SetData(3.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 84)

        # now change a dependent field
        control_field.GetContainerExpressions()[0].SetData(3.0)
        value = self.response_routine.CalculateValue(control_field)
        self.assertEqual(value, 66)

    def test_CalculateGradient(self):
        control_field = self.master_control.GetEmptyControlFields()
        control_field.Read(Kratos.DENSITY)
        value = self.response_routine.CalculateValue(control_field)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()