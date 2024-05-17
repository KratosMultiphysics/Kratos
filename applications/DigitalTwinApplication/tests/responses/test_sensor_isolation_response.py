import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.responses.sensor_isolation_response import SensorIsolationResponse
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorIsolationResponse(UnitTest.TestCase):
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

        cls.sensor_model_part.GetNode(1).SetValue(KratosDT.SENSOR_STATUS, 1)
        cls.sensor_model_part.GetNode(2).SetValue(KratosDT.SENSOR_STATUS, 0.8)
        cls.sensor_model_part.GetNode(3).SetValue(KratosDT.SENSOR_STATUS, 0.7)
        cls.sensor_model_part.GetNode(4).SetValue(KratosDT.SENSOR_STATUS, 0.6)

        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ],
            "isolation_distance": 10.0,
            "beta"              : 300.0
        }""")
        cls.response = SensorIsolationResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.response.CalculateValue(), 0.525)

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
            print(fd_sensitivity, analytical_gradient[i])
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i])

if __name__ == '__main__':
    UnitTest.main()