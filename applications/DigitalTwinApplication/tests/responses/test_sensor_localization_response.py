import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.responses.sensor_localization_response import SensorLocalizationResponse
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorLocalizationResponse(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        """
          (0,0)
            1------2------3------4
            |      |      |      |
            |   1  |  2   |  3   |
            |      |      |      |
            5------6------7------8
            |      |      |      |
            |   4  |  5   |  6   |
            |      |      |      |
            9-----10-----11-----12
                                (6,4)
        """

        cls.model = Kratos.Model()
        cls.mask_model_part = cls.model.CreateModelPart("mask")
        cls.mask_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.mask_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.mask_model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.mask_model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.mask_model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.mask_model_part.CreateNewProperties(1)

        cls.mask_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop)
        cls.mask_model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop)

        for node in cls.mask_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [1, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [3, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [3, 3, 0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [5, 3, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.optimization_problem = OptimizationProblem()

        cls.sensors = GetSensors(cls.mask_model_part, parameters)
        ComponentDataView("sensor", cls.optimization_problem).GetUnBufferedData().SetValue("list_of_sensors", cls.sensors)
        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        for i, sensor in enumerate(cls.sensors):
            loc = sensor.GetLocation()
            node: Kratos.Node = cls.sensor_model_part.CreateNewNode(i + 1, loc[0], loc[1], loc[2])
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 2) / 2)

        elem_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)

        elem_np = np.array([1, 1, 0, 1, 1, 0], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[0].AddElementExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([1, 1, 0, 1, 0, 1], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[1].AddElementExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([0, 0, 0, 1, 0, 0], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[2].AddElementExpression("mask_exp", elem_exp.Clone())

        elem_np = np.array([0, 0, 1, 0, 1, 1], dtype=np.float64)
        Kratos.Expression.CArrayExpressionIO.Read(elem_exp, elem_np)
        cls.sensors[3].AddElementExpression("mask_exp", elem_exp.Clone())

        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ],
            "beta"       : 500
        }""")
        cls.sensor_mask_status = KratosDT.MaskUtils.SensorElementMaskStatus(cls.sensor_model_part, [sensor.GetElementExpression("mask_exp") for sensor in cls.sensors])
        cls.sensor_mask_status_kd_tree = KratosDT.MaskUtils.SensorElementMaskStatusKDTree(cls.sensor_mask_status, 4)
        ComponentDataView("sensor", cls.optimization_problem).GetUnBufferedData().SetValue("sensor_mask_status_kd_tree", cls.sensor_mask_status_kd_tree)
        cls.response = SensorLocalizationResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

        cls.sensor_mask_status.Update()
        cls.sensor_mask_status_kd_tree.Update()

    def test_CalculateValue1(self):
        self.sensor_model_part.GetNode(1).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosDT.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosDT.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1 and m3

        c_1 = e1, e2, e5
        c_2 = e1, e2, e5
        c_3 = e3, e6
        c_4 = e4
        c_5 = e1, e2, e5
        c_6 = e3, e6
        """

        cluster_element_ids_list = [
            [1, 2, 5],
            [1, 2, 5],
            [3, 6],
            [4],
            [1, 2, 5],
            [3, 6]
        ]
        def __calculate_response_value(cluster_element_ids_list: 'list[list[int]]') -> float:
            total_domain_size_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
            Kratos.Expression.DomainSizeExpressionIO.Read(total_domain_size_exp)
            total_domain_size = Kratos.Expression.Utils.Sum(total_domain_size_exp)

            response_value = 0.0
            for cluster_element_ids in cluster_element_ids_list:
                cluster_size = 0.0
                for element_id in cluster_element_ids:
                    cluster_size += self.mask_model_part.GetElement(element_id).GetGeometry().DomainSize()
                cluster_size /= total_domain_size
                response_value += np.exp(cluster_size * 4)
            return np.log(response_value) / 4.0

        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()
        self.assertAlmostEqual(0.5, self.response.CalculateValue())

        self.sensor_model_part.GetNode(1).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(3).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosDT.SENSOR_STATUS, 1)

        """
        status  1   1   1   1
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1, m2, m3, m4

        c_1 = e1, e2
        c_2 = e1, e2
        c_3 = e3
        c_4 = e4
        c_5 = e5
        c_6 = e6
        """

        cluster_element_ids_list = [
            [1, 2],
            [1, 2],
            [3],
            [4],
            [5],
            [6]
        ]

        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()
        self.assertAlmostEqual(1 / 3, self.response.CalculateValue())

    def test_CalculateValue2(self):
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)
        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()
        self.assertAlmostEqual(self.response.CalculateValue(), 0.5)

    def test_CalculateValue3(self):
        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ],
            "beta": 300
        }""")
        response = SensorLocalizationResponse("test1", self.model, params, self.optimization_problem)
        response.Initialize()

        self.sensor_model_part.GetNode(1).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosDT.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosDT.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosDT.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1

        So the clustering will be based on m1 and m3

        c_1 = e1, e2, e5
        c_2 = e1, e2, e5
        c_3 = e3, e6
        c_4 = e4
        c_5 = e1, e2, e5
        c_6 = e3, e6
        """
        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()
        # the exact max cluster size ratio is 0.5
        self.assertAlmostEqual(response.CalculateValue(), 0.5)


    def test_CalculateGradient2(self):
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 4) / 4)
        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()
        ref_value = self.response.CalculateValue()
        collective_exp = KratosOA.CollectiveExpression()
        collective_exp.Add(Kratos.Expression.NodalExpression(self.sensor_model_part))
        self.response.CalculateGradient({KratosDT.SENSOR_STATUS: collective_exp})
        analytical_gradient = collective_exp.GetContainerExpressions()[0].Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) + delta)
            self.sensor_mask_status.Update()
            self.sensor_mask_status_kd_tree.Update()
            fd_sensitivity = (self.response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i], 5)

if __name__ == '__main__':
    UnitTest.main()