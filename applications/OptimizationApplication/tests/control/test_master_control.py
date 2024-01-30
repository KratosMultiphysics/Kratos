import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression

class TestMassterControl(kratos_unittest.TestCase):
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
            "control_variable_name" : "DENSITY"
        }""")
        cls.properties_control_2 = MaterialPropertiesControl("control2", cls.model, parameters)

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test3"],
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

    def test_GetListOfControls(self):
        self.assertEqual([self.properties_control_1, self.properties_control_2, self.properties_control_3, self.properties_control_4], self.master_control.GetListOfControls())

    def test_GetPhysicalKratosVariableCollectiveExpressionsMap(self):
        result = self.master_control.GetPhysicalKratosVariableCollectiveExpressionsMap()
        self.assertEqual([Kratos.DENSITY, Kratos.THICKNESS], list(result.keys()))

        density_collective_expression = result[Kratos.DENSITY]
        density_container_expression_model_part_names = []
        for container_expression in density_collective_expression.GetContainerExpressions():
            self.assertTrue(isinstance(container_expression, Kratos.Expression.ElementExpression))
            density_container_expression_model_part_names.append(container_expression.GetModelPart().FullName())

        self.assertEqual(
            ["test1", "test2", "test3"],
            density_container_expression_model_part_names)

        thickness_collective_expression = result[Kratos.THICKNESS]
        thickness_container_expression_model_part_names = []
        for container_expression in thickness_collective_expression.GetContainerExpressions():
            self.assertTrue(isinstance(container_expression, Kratos.Expression.ElementExpression))
            thickness_container_expression_model_part_names.append(container_expression.GetModelPart().FullName())

        self.assertEqual(
            ["test3"],
            thickness_container_expression_model_part_names)

    def test_GetEmptyField(self):
        empty_control_fields = self.master_control.GetEmptyField()
        container_expression_model_part_names = []
        for container_expression in empty_control_fields.GetContainerExpressions():
            self.assertTrue(isinstance(container_expression, Kratos.Expression.ElementExpression))
            self.assertEqual(Kratos.Expression.Utils.NormInf(container_expression), 0.0)
            container_expression_model_part_names.append(container_expression.GetModelPart().FullName())

        self.assertEqual(
            ["test1", "test2", "test3", "test3"],
            container_expression_model_part_names)

    def test_MapGradient(self):
        result = self.master_control.GetPhysicalKratosVariableCollectiveExpressionsMap()
        mapped_gradients = self.master_control.MapGradient(result)

        for i, control in enumerate(self.master_control.GetListOfControls()):
            self.assertTrue(IsSameContainerExpression(mapped_gradients.GetContainerExpressions()[i], control.GetEmptyField()))

        # now check without any gradients so the master control should fill them with zeros.
        mapped_gradients = self.master_control.MapGradient({})
        for i, control in enumerate(self.master_control.GetListOfControls()):
            mapped_expression = mapped_gradients.GetContainerExpressions()[i]
            self.assertTrue(IsSameContainerExpression(mapped_expression, control.GetEmptyField()))
            self.assertEqual(Kratos.Expression.Utils.NormInf(mapped_expression), 0.0)

    def test_Update(self):
        update = self.master_control.GetEmptyField()

        # assigning density for all the mapped gradients
        KratosOA.CollectiveExpressionIO.Read(update, KratosOA.CollectiveExpressionIO.PropertiesVariable(Kratos.DENSITY))

        # now updating density should not do anything for density controls, but thickness should be updated
        # checking for that
        updated_status = self.master_control.Update(update)
        for k, v in updated_status.items():
            if k.GetPhysicalKratosVariables() == [Kratos.THICKNESS]:
                self.assertTrue(v)
            else:
                self.assertFalse(v)

        update *= 2
        # now everyhing should be updated
        updated_status = self.master_control.Update(update)
        self.assertTrue(all(updated_status.values()))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()