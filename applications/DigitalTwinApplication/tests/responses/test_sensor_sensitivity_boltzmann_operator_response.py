import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.responses.sensor_sensitivity_boltzmann_operator_response import SensorSensitivityBoltzmannOperatorResponse


class TestSensorDistanceBoltzmannOperatorResponseUtils(UnitTest.TestCase):
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

        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        cls.sensor_model_part.CreateNewNode(1, 0, 0, 0)
        cls.sensor_model_part.CreateNewNode(2, 0, 0, 0)
        cls.sensor_model_part.CreateNewNode(3, 0, 0, 0)
        cls.sensor_model_part.CreateNewNode(4, 0, 0, 0)

        """

        theta   0.5     1.0     0.0     0.5
        m       m1      m2      m3      m4      sum
        e1      0.1     0.2     0.3     0.9     (0.05+0.2+0.45)/4=0.175
        e2      0.2     0.0     0.5     0.1     (0.1+0.05)/4     =0.0375
        e3      0.3     0.0     0.7     0.1     (0.15+0.05)/4    =0.05
        e4      0.9     0.2     0.9     0.1     (0.45+0.2+0.05)/4=0.175
        e5      0.5     0.1     0.1     0.1     (0.25+0.1+0.05)/4=0.1
        e6      0.1     0.2     0.0     0.2     (0.05+0.5)/4     =0.0375

        """

        cls.list_of_sensor_sensitivity_distributions = []
        mask = np.array([0.1, 0.2, 0.3, 0.9, 0.5, 0.1], dtype=np.float64)
        mask_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
        Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
        cls.list_of_sensor_sensitivity_distributions.append(mask_exp.Clone())

        mask = np.array([0.2, 0.0, 0.0, 0.2, 0.1, 0.2], dtype=np.float64)
        mask_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
        Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
        cls.list_of_sensor_sensitivity_distributions.append(mask_exp.Clone())

        mask = np.array([0.3, 0.5, 0.7, 0.9, 0.1, 0.0], dtype=np.float64)
        mask_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
        Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
        cls.list_of_sensor_sensitivity_distributions.append(mask_exp.Clone())

        mask = np.array([0.9, 0.1, 0.1, 0.1, 0.1, 0.2], dtype=np.float64)
        mask_exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
        Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
        cls.list_of_sensor_sensitivity_distributions.append(mask_exp.Clone())

    def test_CalculateValue1(self):
        # checking for maximum
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)

        response_utils = KratosDT.SensorElementSensitivityBoltzmannOperatorResponseUtils(self.sensor_model_part, self.list_of_sensor_sensitivity_distributions, 400)
        self.assertAlmostEqual(response_utils.CalculateValue(), 0.175)

    def test_CalculateValue2(self):
        # checking for minimum
        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)

        response_utils = KratosDT.SensorElementSensitivityBoltzmannOperatorResponseUtils(self.sensor_model_part, self.list_of_sensor_sensitivity_distributions, -900)
        self.assertAlmostEqual(response_utils.CalculateValue(), 0.0375, 6)

    def test_CalculateGradient1(self):
        response_utils = KratosDT.SensorElementSensitivityBoltzmannOperatorResponseUtils(self.sensor_model_part, self.list_of_sensor_sensitivity_distributions, -3)

        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosDT.SENSOR_STATUS, (node.Id % 3) / 2)

        ref_value = response_utils.CalculateValue()
        analytical_gradient = response_utils.CalculateGradient().Evaluate()

        delta = 1e-9
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) + delta)
            fd_gradient = (response_utils.CalculateValue() - ref_value) / delta
            node.SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_gradient, analytical_gradient[i], 5)

class TestSensorSensitivityBoltzmannOperatorResponse(UnitTest.TestCase):
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
            exp = Kratos.Expression.ElementExpression(cls.mask_model_part)
            Kratos.Expression.CArrayExpressionIO.Read(exp, np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]) * (i + 1))
            sensor.AddElementExpression("sensitivity", exp)

        params = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "sensors"
            ],
            "sensitivity_expression_name": "sensitivity",
            "beta": -300.0
        }""")
        cls.response = SensorSensitivityBoltzmannOperatorResponse("test", cls.model, params, cls.optimization_problem)
        cls.response.Initialize()

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.response.CalculateValue(), 0.1125)

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