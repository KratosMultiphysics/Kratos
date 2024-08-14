import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSensorMaskStatus(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.sensor_model_part = cls.model.CreateModelPart("senors")
        cls.mask_model_part = cls.model.CreateModelPart("masks")

        number_of_entities = 100
        for i in range(number_of_entities):
            cls.mask_model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)

        number_of_sensors = 50
        cls.sensor_masks_list: 'list[Kratos.Expression.NodalExpression]' = []
        for i in range(number_of_sensors):
            node: Kratos.Node = cls.sensor_model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)
            node.SetValue(KratosSI.SENSOR_STATUS, (i + 15) % 15 / 14)
            for mask_node in cls.mask_model_part.Nodes:
                mask_node.SetValue(Kratos.PRESSURE, (i + mask_node.Id) % 15 / 14)

            nodal_exp = Kratos.Expression.NodalExpression(cls.mask_model_part)
            Kratos.Expression.VariableExpressionIO.Read(nodal_exp, Kratos.PRESSURE, False)
            cls.sensor_masks_list.append(nodal_exp)

        cls.sensor_mask_status = KratosSI.SensorMaskStatus(cls.sensor_model_part,   cls.sensor_masks_list)
        cls.sensor_mask_status.Update()

    def test_GetMasks(self):
        masks = self.sensor_mask_status.GetMasks()

        self.assertEqual(masks.Size1(), self.mask_model_part.NumberOfNodes())
        self.assertEqual(masks.Size2(), self.sensor_model_part.NumberOfNodes())

        for j, _ in enumerate(self.sensor_model_part.Nodes):
            np_exp = self.sensor_masks_list[j].Evaluate()
            for i, _ in enumerate(self.mask_model_part.Nodes):
                self.assertEqual(masks[i, j], np_exp[i])

    def test_GetMaskStatuses(self):
        mask_statuses =  self.sensor_mask_status.GetMaskStatuses()

        self.assertEqual(mask_statuses.Size1(), self.mask_model_part.NumberOfNodes())
        self.assertEqual(mask_statuses.Size2(), self.sensor_model_part.NumberOfNodes())

        for j, sensor_node_j in enumerate(self.sensor_model_part.Nodes):
            np_exp = self.sensor_masks_list[j].Evaluate()
            for i, _ in enumerate(self.mask_model_part.Nodes):
                self.assertEqual(mask_statuses[i, j], np_exp[i] * sensor_node_j.GetValue(KratosSI.SENSOR_STATUS))

if __name__ == '__main__':
    UnitTest.main()