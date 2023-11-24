from math import sqrt
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors


class TestAdjointDisplacementSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        params = Kratos.Parameters("""{
            "model_part_name"        : "Test",
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": false
        }""")
        #cls.adjoint_sensor = KratosDT.Sensors.AdjointDisplacementSensor(cls.model, params)

    def test_CalculateGradient(self):
        residual_sensor = KratosDT.Sensors.MeasurementResidualPNormResponseFunction(2.0)

        params_list = Kratos.Parameters("""{
            "list_of_sensors": [
                {
                        "type"         : "displacement_sensor",
                        "name"         : "s1",
                        "value"        : 0,
                        "location"     : [0.3333, 0.3333, 0.0],
                        "direction"    : [1.0, 0.0, 0.0],
                        "weight"       : 1.0,
                        "variable_data": {}
                },
                {
                        "type"         : "displacement_sensor",
                        "name"         : "s2",
                        "value"        : 0,
                        "location"     : [0.666, 0.666, 0.0],
                        "direction"    : [1.0, 0.0, 0.0],
                        "weight"       : 1.0,
                        "variable_data": {}
                }
            ]
        }""")

        sensors_list = GetSensors(self.model_part, params_list["list_of_sensors"].values())
        for sensor in sensors_list:
            sensor.SetValue(KratosDT.SENSOR_MEASURED_VALUE, 1.0)
            residual_sensor.AddSensor(sensor)

        ref_value = residual_sensor.CalculateValue(self.model_part)
        print("ref_value", ref_value)
        delta = 1e-6

        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        element = self.model_part.GetElement(1)
        residual_sensor.CalculateValue(self.model_part)
        residual_sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
        print(response_sensitivities)
        for i, node in enumerate(element.GetGeometry()):
            print(i)
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
            perturbed_value = residual_sensor.CalculateValue(self.model_part)
            print("perturbed_value", perturbed_value)
            sensitivity = (perturbed_value - ref_value) / delta
            print(sensitivity, response_sensitivities[i * 6])
            # self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6], 6)
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
            perturbed_value = residual_sensor.CalculateValue(self.model_part)
            sensitivity = (perturbed_value - ref_value) / delta
            # print(sensitivity, response_sensitivities[i * 6 ] + 1)
            self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1], 6)
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


if __name__ == '__main__':
    UnitTest.main()