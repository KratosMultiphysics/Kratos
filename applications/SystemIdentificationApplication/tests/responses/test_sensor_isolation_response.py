import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.sensor_isolation_response import SensorIsolationResponse

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

        cls.sensors = GetSensors(cls.mask_model_part, parameters)
        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        for i, sensor in enumerate(cls.sensors):
            loc = sensor.GetLocation()
            node: Kratos.Node = cls.sensor_model_part.CreateNewNode(i + 1, loc[0], loc[1], loc[2])
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 3) / 2)

        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ],
            "isolation_radius": 10.0
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.response = SensorIsolationResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue(self):
        distance = 0.0
        for node_i in self.sensor_model_part.Nodes:
            for node_j in self.sensor_model_part.Nodes:
                if node_i.Id != node_j.Id:
                    current_distance = ((node_i.X - node_j.X) ** 2 + (node_i.Y - node_j.Y) ** 2 + (node_i.Z - node_j.Z) ** 2) ** (0.5)
                    if current_distance > 10.0:
                        current_distance = 0.0
                    else:
                        x = 1.0 - current_distance / 10.0
                        current_distance = x*x*(3.0-2.0*x)
                    distance += node_i.GetValue(KratosSI.SENSOR_STATUS) * node_j.GetValue(KratosSI.SENSOR_STATUS) * current_distance
        self.assertAlmostEqual(self.response.CalculateValue(), distance * 0.5)

    def test_CalculateGradient(self):
        ref_value = self.response.CalculateValue()
        collective_exp = KratosOA.CollectiveExpression()
        collective_exp.Add(Kratos.Expression.NodalExpression(self.sensor_model_part))
        self.response.CalculateGradient({KratosSI.SENSOR_STATUS: collective_exp})
        analytical_gradient = collective_exp.GetContainerExpressions()[0].Evaluate()

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) + delta)
            fd_sensitivity = (self.response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient[i])

if __name__ == '__main__':
    UnitTest.main()