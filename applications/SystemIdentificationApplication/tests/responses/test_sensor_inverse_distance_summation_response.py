import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.responses.sensor_inverse_distance_summation_response import SensorInverseDistanceSummationResponse

class TestSensorInverseDistanceSummationResponse(UnitTest.TestCase):
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
        cls.domain_model_part = cls.model.CreateModelPart("mask")
        cls.domain_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.domain_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.domain_model_part.CreateNewProperties(1)

        cls.domain_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop)

        for node in cls.domain_model_part.Nodes:
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

        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        cls.sensors = CreateSensors(cls.sensor_model_part, cls.domain_model_part, parameters)

        cls.optimization_problem = OptimizationProblem()
        cls.sensor_group_data = ComponentDataView("sensors", cls.optimization_problem)
        SetSensors(cls.sensor_group_data, cls.sensors)

        for sensor in cls.sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_STATUS, (sensor.GetNode().Id % 3) / 2)

        params = Kratos.Parameters("""{
            "sensor_group_name": "sensors",
            "p_coefficient": 2.0
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.response = SensorInverseDistanceSummationResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue(self):
        distance = 0.0
        for node_i in self.sensor_model_part.Nodes:
            for node_j in self.sensor_model_part.Nodes:
                if node_i.Id != node_j.Id:
                    current_distance = ((node_i.X - node_j.X) ** 2 + (node_i.Y - node_j.Y) ** 2 + (node_i.Z - node_j.Z) ** 2) ** (0.5)
                    distance += ((node_i.GetValue(KratosSI.SENSOR_STATUS) * node_j.GetValue(KratosSI.SENSOR_STATUS)) ** 2) / current_distance
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