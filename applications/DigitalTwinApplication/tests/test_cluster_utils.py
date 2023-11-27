import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from random import shuffle

class TestClusterUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0, 0, 0)
        cls.model_part.CreateNewNode(2, 1, 0, 0)
        cls.model_part.CreateNewNode(3, 1, 1, 0)

        prop = cls.model_part.CreateNewProperties(1)
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        params = Kratos.Parameters("""{
            "list_of_sensors": [
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_1",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_2",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_3",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_4",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_5",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_6",
                    "type": "displacement_sensor"
                },
                {
                    "direction": [1.0, 0.0, 0.0],
                    "location": [0.33, 0.33, 0.33],
                    "name": "disp_x_7",
                    "type": "displacement_sensor"
                }
            ]
        }""")

        cls.sensors = GetSensors(cls.model_part, params["list_of_sensors"].values())

        cls.sensor_views: 'list[KratosDT.Sensors.NodalSensorView]' = []
        # assign dummy expressions and sensor views
        for i, sensor in enumerate(cls.sensors):
            nodal_exp = Kratos.Expression.NodalExpression(cls.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(nodal_exp, numpy.arange(i + 1, i + cls.model_part.NumberOfNodes() + 1, dtype=numpy.float64))
            sensor.AddContainerExpression("test_expression", nodal_exp)
            cls.sensor_views.append(KratosDT.Sensors.NodalSensorView(sensor, "test_expression"))

        cls.distances_matrix: 'list[list[float]]' = []
        for sensor_view_i in cls.sensor_views:
            distances_row: 'list[float]' = []
            for sensor_view_j in cls.sensor_views:
                distances_row.append(KratosOA.ExpressionUtils.InnerProduct(sensor_view_i.GetContainerExpression(), sensor_view_j.GetContainerExpression()))
            cls.distances_matrix.append(distances_row)

        cls.distances_compressed_matrix: 'list[float]' = []
        for i, row in enumerate(cls.distances_matrix):
            for col in row[i+1:]:
                cls.distances_compressed_matrix.append(col)

        cls.clustering_data = KratosDT.ClusterUtils.NodalSensorViewClusterData(cls.sensor_views, [cls.sensor_views[1], cls.sensor_views[3], cls.sensor_views[5]])
        cls.clustering_data.AddDistances("dist", cls.distances_compressed_matrix)

    def test_ClusterIds(self):
        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        self.assertEqual(cluster.Id, 1)

        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(2, self.clustering_data)
        self.assertEqual(cluster.Id, 2)

    def test_GetClusterSensorViews(self):
        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        shuffeled_list = list(self.sensor_views)
        shuffle(shuffeled_list)
        cluster.SetSensorViews(shuffeled_list)
        self.assertEqual(cluster.GetSensorViews(), self.sensor_views)

        sub_cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        sub_cluster.SetSensorViews([self.sensor_views[5], self.sensor_views[3], self.sensor_views[2], self.sensor_views[1]])
        self.assertEqual(sub_cluster.GetSensorViews(), [self.sensor_views[1], self.sensor_views[2], self.sensor_views[3], self.sensor_views[5]])

    def test_GetClusterDistances(self):
        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        cluster.SetSensorViews(self.sensor_views)
        distances = cluster.GetDistances("dist")
        self.assertEqual(distances, self.distances_compressed_matrix)


        sub_cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        sub_cluster.SetSensorViews([self.sensor_views[5], self.sensor_views[3], self.sensor_views[2], self.sensor_views[1]])
        distances = sub_cluster.GetDistances("dist")
        sub_indices = [1, 2, 3, 5]
        for i, sub_index_i in enumerate(sub_indices):
            for temp_j, sub_index_j in enumerate(sub_indices[i+1:]):
                j = temp_j + i + 1
                sub_distance = distances[len(sub_indices) * i + j - ((i+2)*(i+1)//2)]
                self.assertEqual(self.distances_matrix[sub_index_i][sub_index_j], sub_distance)

    def test_GetEntities(self):
        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        cluster.SetSensorViews(self.sensor_views)
        cluster_entities = cluster.GetEntities()
        self.assertEqual(len(self.model_part.Nodes), len(cluster_entities))
        for ref_entity, cluster_entity  in zip(self.model_part.Nodes, cluster_entities):
            self.assertEqual(ref_entity, cluster_entity)

        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        cluster.SetSensorViews([self.sensor_views[0], self.sensor_views[6]])
        cluster_entities = cluster.GetEntities()
        self.assertEqual(len(cluster_entities), 0)

        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        cluster.SetSensorViews([self.sensor_views[3], self.sensor_views[6], self.sensor_views[0], self.sensor_views[1]])
        cluster_entities = cluster.GetEntities()
        self.assertEqual(len(cluster_entities), 2)
        for ref_entity, cluster_entity  in zip([self.model_part.GetNode(1), self.model_part.GetNode(2)], cluster_entities):
            self.assertEqual(ref_entity, cluster_entity)

        cluster = KratosDT.ClusterUtils.NodalSensorViewCluster(1, self.clustering_data)
        cluster.SetSensorViews([self.sensor_views[5], self.sensor_views[6], self.sensor_views[0], self.sensor_views[1]])
        cluster_entities = cluster.GetEntities()
        self.assertEqual(len(cluster_entities), 2)
        for ref_entity, cluster_entity  in zip([self.model_part.GetNode(1), self.model_part.GetNode(3)], cluster_entities):
            self.assertEqual(ref_entity, cluster_entity)

if __name__ == '__main__':
    UnitTest.main()