import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.responses.sensor_distance_variance_response import SensorDistanceVarianceResponse
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView


class TestSensorDistanceVarianceUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.CreateNewNode(1, 0, 0, 0)
        cls.model_part.CreateNewNode(2, 3, 0, 0)
        cls.model_part.CreateNewNode(3, 5, 0, 0)
        cls.model_part.CreateNewNode(4, 5, 5, 0)
        cls.model_part.CreateNewNode(5, 10, 20, 0)
        cls.model_part.CreateNewNode(6, 10, 40, 0)

        cls.response_utils = KratosDT.SensorDistanceVarianceUtils(cls.model_part)
        cls.response_utils.Initialize()

    def test_CalculateValue1(self):
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, (i % 3) / 2)

        n = self.model_part.NumberOfNodes()
        distances = []
        for i in range(n):
            n_i: Kratos.Node = self.model_part.GetNode(i + 1)
            for j in range(i + 1, n):
                n_j: Kratos.Node = self.model_part.GetNode(j + 1)
                distance = ((n_i.X - n_j.X) ** 2 + (n_i.Y - n_j.Y) ** 2 + (n_i.Z - n_j.Z) ** 2) ** 0.5
                distances.append(distance * n_i.GetValue(KratosDT.SENSOR_STATUS) * n_j.GetValue(KratosDT.SENSOR_STATUS))

        self.assertAlmostEqual(self.response_utils.CalculateValue(), np.var(distances))

    def test_CalculateGradient(self):
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, (i % 3) / 2)

        ref_value = self.response_utils.CalculateValue()
        analytical_gradient = self.response_utils.CalculateGradient().Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) + delta)
            fd_gradient = (self.response_utils.CalculateValue() - ref_value) / delta
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_gradient, analytical_gradient[i], 5)

class TestSensorDistanceVarianceResponse(UnitTest.TestCase):
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
                "location"     : [1.0001, 1, 0.0],
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
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)

        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ]
        }""")
        cls.response = SensorDistanceVarianceResponse("test", cls.model, params)
        cls.response.Initialize()

    def test_CalculateValue(self):
        n = self.sensor_model_part.NumberOfNodes()
        distances = []
        for i in range(n):
            n_i: Kratos.Node = self.sensor_model_part.GetNode(i + 1)
            for j in range(i + 1, n):
                n_j: Kratos.Node = self.sensor_model_part.GetNode(j + 1)
                distance = ((n_i.X - n_j.X) ** 2 + (n_i.Y - n_j.Y) ** 2 + (n_i.Z - n_j.Z) ** 2) ** 0.5
                distances.append(distance * n_i.GetValue(KratosDT.SENSOR_STATUS) * n_j.GetValue(KratosDT.SENSOR_STATUS))


        self.assertAlmostEqual(self.response.CalculateValue(), np.var(distances))

    def test_CalculateGradient(self):
        ref_value = self.response.CalculateValue()
        collective_exp = KratosOA.CollectiveExpression()
        collective_exp.Add(Kratos.Expression.NodalExpression(self.sensor_model_part))
        self.response.CalculateGradient({KratosDT.SENSOR_STATUS: collective_exp})
        analytical_gradient = collective_exp.GetContainerExpressions()[0].Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) + delta)
            fd_sensitivity = (self.response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i], 6)

if __name__ == '__main__':
    UnitTest.main()