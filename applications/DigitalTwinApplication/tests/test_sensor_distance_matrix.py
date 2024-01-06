import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSensorUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0, 0, 0)
        cls.model_part.CreateNewNode(2, 100, 0, 0)
        cls.model_part.CreateNewNode(3, 0, 100, 0)

        prop = cls.model_part.CreateNewProperties(1)
        cls.element = cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        cls.n = 10
        cls.sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
        for i in range(cls.n):
            sensor = KratosDT.Sensors.DisplacementSensor(f"{i}", Kratos.Point(i, i+1, i+2), [1, 0, 0], cls.element, 1.0)
            cls.sensors_list.append(sensor)

        cls.sensor_distance_matrix = KratosDT.SensorDistanceMatrix(cls.sensors_list)

    def test_GetDistance1(self):
        for i in range(self.n):
            for j in range(i, self.n):
                distance = self.sensors_list[i].GetLocation() - self.sensors_list[j].GetLocation()
                distance = (distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2) ** 0.5
                self.assertAlmostEqual(distance, self.sensor_distance_matrix.GetDistance(self.sensors_list[i], self.sensors_list[j]))

    def test_GetDistance2(self):
        sensor_list_1 = [self.sensors_list[i] for i in [1,2,3,4]]
        sensor_list_2 = [self.sensors_list[i] for i in [0,7,5,9,8]]
        distance_matrix = self.sensor_distance_matrix.GetDistance(sensor_list_1, sensor_list_2)

        for i, sensor_1 in enumerate(sensor_list_1):
            for j, sensor_2 in enumerate(sensor_list_2):
                distance = sensor_1.GetLocation() - sensor_2.GetLocation()
                distance = (distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2) ** 0.5
                self.assertAlmostEqual(distance, distance_matrix[i, j])

    def test_GetDistance3(self):
        sensor_list_1 = [self.sensors_list[i] for i in [1,2,3,4]]
        sensor_list_2 = [self.sensors_list[i] for i in [0,7,5,9,8]]
        sensor_list_2.append(KratosDT.Sensors.DisplacementSensor("dummy", Kratos.Point(1, 1, 1), [1, 0, 0], self.element, 1.0))
        distance_matrix = self.sensor_distance_matrix.GetDistance(sensor_list_1, sensor_list_2)

        for i, sensor_1 in enumerate(sensor_list_1):
            for j, sensor_2 in enumerate(sensor_list_2[:-1]):
                distance = sensor_1.GetLocation() - sensor_2.GetLocation()
                distance = (distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2) ** 0.5
                self.assertAlmostEqual(distance, distance_matrix[i, j])

            self.assertEqual(distance_matrix[i, len(sensor_list_2) - 1], -1)

    def test_MaxDistances(self):
        m1 = Kratos.Matrix(10, 10, 0.0)
        m2 = Kratos.Matrix(10, 10, 0.0)

        for i in range(10):
            for j in range(10):
                m1[i, j] = (i + 1) * (j + 1) * (-1) ** (i + j)
                m2[i, j] = (i + 1) * (j + 1) * (-1) ** (i + j + 1)

        m3 = KratosDT.SensorDistanceMatrix.MaxDistances(m1, m2)
        for i in range(10):
            for j in range(10):
                self.assertEqual(m3[i, j], max(m1[i, j], m2[i, j]))

if __name__ == '__main__':
    UnitTest.main()
