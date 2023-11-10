from math import sqrt
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.KratosUnittest as UnitTest


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
        cls.adjoint_sensor = KratosDT.Sensors.AdjointDisplacementSensor(cls.model, params)

    def test_CalculateValue(self):
        specification = KratosDT.Sensors.SensorSpecification("test", 1)
        specification.SetLocation([1/3, 1/3, 0.0])
        specification.SetValue(KratosDT.SENSOR_WEIGHT, 1.0)
        specification.SetValue(KratosDT.SENSOR_ELEMENT_ID, 1)
        specification.SetValue(KratosDT.SENSOR_DIRECTION, [1, 0, 0])
        self.adjoint_sensor.SetSensorSpecification(specification)
        value = self.adjoint_sensor.CalculateValue(self.model_part)
        self.assertAlmostEqual(value, 7/3)

        specification.SetValue(KratosDT.SENSOR_DIRECTION, [1/sqrt(2.0), 1/sqrt(2.0), 0])
        self.adjoint_sensor.SetSensorSpecification(specification)
        value = self.adjoint_sensor.CalculateValue(self.model_part)
        self.assertAlmostEqual(value, 4.006938426723769)

        specification.SetLocation([2/3, 2/3, 0])
        specification.SetValue(KratosDT.SENSOR_ELEMENT_ID, 2)
        specification.SetValue(KratosDT.SENSOR_DIRECTION, [1, 0, 0])
        self.adjoint_sensor.SetSensorSpecification(specification)
        value = self.adjoint_sensor.CalculateValue(self.model_part)
        self.assertAlmostEqual(value, 3.0)

        specification.SetValue(KratosDT.SENSOR_DIRECTION, [1/sqrt(2.0), 1/sqrt(2.0), 0])
        self.adjoint_sensor.SetSensorSpecification(specification)
        value = self.adjoint_sensor.CalculateValue(self.model_part)
        self.assertAlmostEqual(value, 4.949747468305832)

    def test_CalculateGradient(self):
        specification = KratosDT.Sensors.SensorSpecification("test", 1)
        specification.SetLocation([1/3, 1/3, 0.0])
        specification.SetValue(KratosDT.SENSOR_WEIGHT, 1.0)
        specification.SetValue(KratosDT.SENSOR_ELEMENT_ID, 1)
        specification.SetValue(KratosDT.SENSOR_DIRECTION, [1/sqrt(2.0), 1/sqrt(2.0), 0])
        self.adjoint_sensor.SetSensorSpecification(specification)

        ref_value = self.adjoint_sensor.CalculateValue(self.model_part)
        delta = 1e-8

        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        element = self.model_part.GetElement(1)
        self.adjoint_sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
        for i, node in enumerate(element.GetGeometry()):
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
            perturbed_value = self.adjoint_sensor.CalculateValue(self.model_part)
            sensitivity = (perturbed_value - ref_value) / delta
            self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6])
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
            perturbed_value = self.adjoint_sensor.CalculateValue(self.model_part)
            sensitivity = (perturbed_value - ref_value) / delta
            self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1])
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


if __name__ == '__main__':
    UnitTest.main()