import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSubDistances

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

    def test_GetSubDistances(self):
        distances_matrix = [
            [0,  2,  3,  4,  5,  6],
            [2,  0,  7,  8,  9, 10],
            [3,  7,  0, 11, 12, 13],
            [4,  8, 11,  0, 14, 15],
            [5,  9, 12, 14,  0, 16],
            [6, 10, 13, 15, 16,  0]
        ]


        def get_distances_vector(matrix_of_distances):
            distances_vector: 'list[int]' = []
            for i, row in enumerate(matrix_of_distances):
                for col in row[i+1:]:
                    distances_vector.append(col)
            return distances_vector

        sub_distances = GetSubDistances(get_distances_vector(distances_matrix), [1,3,4,5])
        sub_matrix = []
        for i in [1,3,4,5]:
            row = []
            for j in [1,3,4,5]:
                row.append(distances_matrix[i][j])
            sub_matrix.append(row)

        self.assertEqual(sub_distances, get_distances_vector(sub_matrix))

    def test_AssignConsecutiveSensorIds(self):
        n = 10
        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
        for i in range(n):
            sensor = KratosDT.Sensors.DisplacementSensor(f"{i}", Kratos.Point(i, i+1, i+2), [1, 0, 0], self.element, 1.0)
            sensors_list.append(sensor)

        KratosDT.SensorUtils.AssignConsecutiveSensorIds(sensors_list, KratosDT.SENSOR_ID)

        for i, sensor in enumerate(sensors_list):
            self.assertEqual(i+1, sensor.GetValue(KratosDT.SENSOR_ID))

    def test_AssignSensorNeighbours(self):
        n = 10
        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
        for i in range(n):
            sensor = KratosDT.Sensors.DisplacementSensor(f"{i}", Kratos.Point(i, i+1, i+2), [1, 0, 0], self.element, 1.0)
            sensors_list.append(sensor)

        KratosDT.SensorUtils.AssignConsecutiveSensorIds(sensors_list, KratosDT.SENSOR_ID)
        KratosDT.SensorUtils.AssignSensorNeighbours(sensors_list, 2, 100, Kratos.INITIAL_STRAIN, KratosDT.SENSOR_ID)

        for sensor in sensors_list:
            neighbour_ids: 'list[int]' = []
            for other_sensor in sensors_list:
                distance = sensor.GetLocation() - other_sensor.GetLocation()
                distance = distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2
                if distance < 4.0:
                    neighbour_ids.append(other_sensor.GetValue(KratosDT.SENSOR_ID))

            self.assertVectorAlmostEqual(neighbour_ids, sensor.GetValue(Kratos.INITIAL_STRAIN))

    def test_ComputeSensorRobustness(self):
        n = 5
        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
        expressions_list: 'list[Kratos.Expression.NodalExpression]' = []
        for i in range(n):
            sensor = KratosDT.Sensors.DisplacementSensor(f"{i}", Kratos.Point(i, i+1, i+2), [1, 0, 0], self.element, 1.0)
            np_expression = numpy.arange(i, i+3, dtype=numpy.float64)

            expression = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(expression, np_expression)

            sensors_list.append(sensor)
            expressions_list.append(expression)

        KratosDT.SensorUtils.AssignConsecutiveSensorIds(sensors_list, KratosDT.SENSOR_ID)
        KratosDT.SensorUtils.AssignSensorNeighbours(sensors_list, 100, 100, Kratos.INITIAL_STRAIN, KratosDT.SENSOR_ID)
        KratosDT.SensorUtils.ComputeSensorRobustness(sensors_list, expressions_list, Kratos.INITIAL_STRAIN, KratosDT.SENSOR_ID, Kratos.PRESSURE)

        for i, sensor in enumerate(sensors_list):
            avg_cosine_distance = 0.0
            i_np_exp = expressions_list[i].Evaluate()
            for j, exp in enumerate(expressions_list):
                if i != j:
                    j_np_exp = exp.Evaluate()
                    avg_cosine_distance += numpy.inner(i_np_exp, j_np_exp) / (numpy.linalg.norm(i_np_exp) * numpy.linalg.norm(j_np_exp))
            self.assertAlmostEqual(avg_cosine_distance / n, sensor.GetValue(Kratos.PRESSURE))

if __name__ == '__main__':
    UnitTest.main()
