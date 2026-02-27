import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl

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

    def test_GetPhysicalKratosVariableCombinedTensorAdaptorsMap(self):
        result = self.master_control.GetPhysicalKratosVariableCombinedTensorAdaptorsMap()
        self.assertEqual([Kratos.DENSITY, Kratos.THICKNESS], list(result.keys()))

        density_cta = result[Kratos.DENSITY]
        density_ta_containers = []
        for ta in density_cta.GetTensorAdaptors():
            self.assertTrue(isinstance(ta.GetContainer(), Kratos.ElementsArray))
            density_ta_containers.append(ta.GetContainer())

        self.assertEqual(
            [self.model["test1"].Elements, self.model["test2"].Elements, self.model["test3"].Elements],
            density_ta_containers)

        thickness_ta = result[Kratos.THICKNESS]
        thickness_containers = []
        for ta in thickness_ta.GetTensorAdaptors():
            self.assertTrue(isinstance(ta.GetContainer(), Kratos.ElementsArray))
            thickness_containers.append(ta.GetContainer())

        self.assertEqual(
            [self.model["test3"].Elements],
            thickness_containers)

    def test_GetEmptyField(self):
        empty_control_fields = self.master_control.GetEmptyField()
        containers = []
        for ta in empty_control_fields.GetTensorAdaptors():
            self.assertTrue(isinstance(ta.GetContainer(), Kratos.ElementsArray))
            self.assertEqual(numpy.linalg.norm(ta.data), 0.0)
            containers.append(ta.GetContainer())

        self.assertEqual(
            [self.model["test1"].Elements, self.model["test2"].Elements, self.model["test3"].Elements, self.model["test3"].Elements],
            containers)

    def test_MapGradient(self):
        result = self.master_control.GetPhysicalKratosVariableCombinedTensorAdaptorsMap()
        mapped_gradients = self.master_control.MapGradient(result)

        for i, control in enumerate(self.master_control.GetListOfControls()):
            self.assertTrue(mapped_gradients.GetTensorAdaptors()[i].GetContainer(), control.GetEmptyField().GetContainer())

        # now check without any gradients so the master control should fill them with zeros.
        mapped_gradients = self.master_control.MapGradient({})
        for i, control in enumerate(self.master_control.GetListOfControls()):
            mapped_ta = mapped_gradients.GetTensorAdaptors()[i]
            self.assertTrue(mapped_ta.GetContainer(), control.GetEmptyField().GetContainer())
            self.assertEqual(numpy.linalg.norm(mapped_ta.data), 0.0)

    def test_Update(self):
        update = self.master_control.GetEmptyField()

        # assigning density for all the mapped gradients
        for ta in update.GetTensorAdaptors():
            KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(ta, Kratos.DENSITY, copy=False).CollectData()
        update.CollectData()

        # now updating density should not do anything for density controls, but thickness should be updated
        # checking for that
        updated_status = self.master_control.Update(update)
        for k, v in updated_status.items():
            if k.GetPhysicalKratosVariables() == [Kratos.THICKNESS]:
                self.assertTrue(v)
            else:
                self.assertFalse(v)

        update.data[:] *= 2.0
        update.StoreData()
        # now everything should be updated
        updated_status = self.master_control.Update(update)
        self.assertTrue(all(updated_status.values()))

if __name__ == "__main__":
    kratos_unittest.main()