from math import fabs
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSensorMaskStatusKDTree(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.sensor_model_part = cls.model.CreateModelPart("senors")
        cls.mask_model_part = cls.model.CreateModelPart("masks")

        number_of_entities = 10
        for i in range(number_of_entities + 1):
            cls.mask_model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)

        properties = cls.mask_model_part.CreateNewProperties(1)
        for i in range(number_of_entities):
            cls.mask_model_part.CreateNewElement("Element3D2N", i + 1,  [i + 1, i + 2], properties)

        number_of_sensors = 5
        cls.sensor_masks_list: 'list[Kratos.Expression.NodalExpression]' = []
        for i in range(number_of_sensors):
            node: Kratos.Node = cls.sensor_model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)
            node.SetValue(KratosSI.SENSOR_STATUS, (i + 15) % 15 / 14)
            for mask_element in cls.mask_model_part.Elements:
                mask_element.SetValue(Kratos.PRESSURE, (i + mask_element.Id) % 15 / 14)

            element_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
            Kratos.Expression.VariableExpressionIO.Read(element_exp, Kratos.PRESSURE)
            cls.sensor_masks_list.append(element_exp)

        cls.sensor_mask_status = KratosSI.SensorMaskStatus(cls.sensor_model_part, cls.sensor_masks_list, 0)
        cls.sensor_mask_status.Update()

        cls.sensor_mask_status_kd_tree = KratosSI.SensorMaskStatusKDTree(cls.sensor_mask_status, 100, 0)
        cls.sensor_mask_status_kd_tree.Update()

    def test_GetSensorMaskStatus(self):
        self.assertEqual(self.sensor_mask_status_kd_tree.GetSensorMaskStatus(), self.sensor_mask_status)

    def test_RadiusSearch(self):
        mask_status = self.sensor_mask_status.GetMaskStatuses()

        radius = 0.03
        for query_point_index in range(mask_status.Size1()):
            ref_neighbours: 'list[tuple[int, float]]' = []
            for i in range(mask_status.Size1()):
                distance = 0.0
                for k in range(mask_status.Size2()):
                    distance += fabs(mask_status[i, k] - mask_status[query_point_index, k])
                if distance <= radius:
                    ref_neighbours.append((i, distance))

            query_point = Kratos.Vector(self.sensor_model_part.NumberOfNodes())
            for i in range(query_point.Size()):
                query_point[i] = mask_status[query_point_index, i]

            neighbours = self.sensor_mask_status_kd_tree.RadiusSearch(query_point, radius)
            neighbours = sorted(neighbours, key=lambda x: x[0])

            self.assertEqual(neighbours, ref_neighbours)

if __name__ == '__main__':
    UnitTest.main()