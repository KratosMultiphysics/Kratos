import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors

class TestDamageDetectionAdjointResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 2.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 3.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 4.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)

        cls.adjoint_response_function = KratosSI.Responses.MeasurementResidualResponseFunction(3.0)

        for i, sensor in enumerate(cls.sensors):
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, i * 15 - 10)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
            cls.adjoint_response_function.AddSensor(sensor)

        cls.adjoint_response_function.Initialize()

    def setUp(self) -> None:
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
        self.ref_value = self.adjoint_response_function.CalculateValue(self.model_part)

    def test_CalculateValue(self):
        value = 0.0
        for sensor in self.sensors:
            value += (0.5 * sensor.GetWeight() * (sensor.CalculateValue(self.model_part) - sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)) ** 2) ** 3.0
        self.assertAlmostEqual(self.ref_value, value ** (1 / 3), 5)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()

        global_fd_x_sensitivities: 'dict[int, float]' = {}
        global_fd_y_sensitivities: 'dict[int, float]' = {}
        analytical_x_sensitivities: 'dict[int, float]' = {}
        analytical_y_sensitivities: 'dict[int, float]' = {}

        delta = 1e-7
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
            perturbed_value = self.adjoint_response_function.CalculateValue(self.model_part)
            global_fd_x_sensitivities[node.Id] = (perturbed_value - self.ref_value) / delta
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
            perturbed_value = self.adjoint_response_function.CalculateValue(self.model_part)
            global_fd_y_sensitivities[node.Id] = (perturbed_value - self.ref_value) / delta
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)

            analytical_x_sensitivities[node.Id] = 0.0
            analytical_y_sensitivities[node.Id] = 0.0

        for element in self.model_part.Elements:
            self.adjoint_response_function.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                analytical_x_sensitivities[node.Id] += response_sensitivities[i * 6]
                analytical_y_sensitivities[node.Id] += response_sensitivities[i * 6 + 1]

        for node_id in global_fd_x_sensitivities.keys():
            self.assertAlmostEqual(global_fd_x_sensitivities[node_id], analytical_x_sensitivities[node_id], 5)
            self.assertAlmostEqual(global_fd_y_sensitivities[node_id], analytical_y_sensitivities[node_id], 5)

class TestDamageDetectionAdjointResponseFunctionErrorThreshold(kratos_unittest.TestCase):
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

        # baseline sensors with a negligible threshold to determine the raw (un-thresholded) sensor errors
        cls.baseline_response, cls.baseline_sensors = cls.__CreateResponseAndSensors(1e-16, "Baseline")
        cls.baseline_response.CalculateValue(cls.model_part)
        cls.raw_errors = [sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR) for sensor in cls.baseline_sensors]

    @classmethod
    def __CreateResponseAndSensors(cls, error_threshold, suffix: str):
        """Creates the 4 test sensors and the response function.

        If error_threshold is None, the "error_threshold" key is omitted from the
        sensor parameters entirely, so DisplacementSensor::GetDefaultParameters()'s
        default value (Sensor::DefaultErrorThreshold) is the one that gets applied.
        """
        sensor_model_part = cls.model.CreateModelPart(f"SensorModelPart{suffix}")

        error_threshold_entry = "" if error_threshold is None else '"error_threshold": %.20g,' % error_threshold

        sensor_settings = [
            ("disp_x_1", "[0.3333333333333, 0.3333333333333, 0.0]", "[1.0, 0.0, 0.0]", 2.0),
            ("disp_x_2", "[0.6666666666667, 0.6666666666667, 0.0]", "[1.0, 0.0, 0.0]", 3.0),
            ("disp_x_3", "[0.3333333333333, 0.3333333333333, 0.0]", "[1.0, 1.0, 0.0]", 4.0),
            ("disp_x_4", "[0.6666666666667, 0.6666666666667, 0.0]", "[1.0, 1.0, 0.0]", 1.0),
        ]

        parameters = [
            Kratos.Parameters("""{
                "type"         : "displacement_sensor",
                "name"         : "%s",
                "value"        : 0,
                "location"     : %s,
                "direction"    : %s,
                "weight"       : %s,
                %s
                "variable_data": {}
            }""" % (name, location, direction, weight, error_threshold_entry))
            for name, location, direction, weight in sensor_settings
        ]

        sensors = CreateSensors(sensor_model_part, cls.model_part, parameters)

        response = KratosSI.Responses.MeasurementResidualResponseFunction(3.0)
        for i, sensor in enumerate(sensors):
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, i * 15 - 10)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
            response.AddSensor(sensor)
        response.Initialize()

        return response, sensors

    def __ExpectedValue(self, sensors, error_threshold: float) -> float:
        value = 0.0
        for sensor in sensors:
            error = sensor.CalculateValue(self.model_part) - sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)
            if abs(error) < error_threshold:
                error = 0.0
            value += (0.5 * sensor.GetWeight() * error ** 2) ** 3.0
        return value ** (1 / 3)

    def test_ThresholdBelowAllErrors(self):
        # threshold lower than every sensor error -> SENSOR_ERROR should keep its original value
        error_threshold = min(abs(e) for e in self.raw_errors) / 2.0

        response, sensors = self.__CreateResponseAndSensors(error_threshold, "BelowAll")
        value = response.CalculateValue(self.model_part)

        self.assertAlmostEqual(value, self.__ExpectedValue(sensors, error_threshold), 8)
        for sensor, raw_error in zip(sensors, self.raw_errors):
            self.assertAlmostEqual(sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR), raw_error, 8)

    def test_ThresholdAboveAllErrors(self):
        # threshold higher than every sensor error -> all SENSOR_ERROR values should be set to 0.0
        error_threshold = max(abs(e) for e in self.raw_errors) + 1.0

        response, sensors = self.__CreateResponseAndSensors(error_threshold, "AboveAll")
        value = response.CalculateValue(self.model_part)

        self.assertAlmostEqual(value, 0.0, 8)
        for sensor in sensors:
            self.assertEqual(sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR), 0.0)

    def test_ThresholdAboveSomeErrors(self):
        # threshold between the two smallest error magnitudes -> some SENSOR_ERROR values are
        # zeroed out while the rest keep their original value
        sorted_abs_errors = sorted(abs(e) for e in self.raw_errors)
        self.assertGreater(sorted_abs_errors[1], sorted_abs_errors[0], "Test setup requires distinct sensor error magnitudes.")
        error_threshold = (sorted_abs_errors[0] + sorted_abs_errors[1]) / 2.0

        response, sensors = self.__CreateResponseAndSensors(error_threshold, "AboveSome")
        value = response.CalculateValue(self.model_part)

        self.assertAlmostEqual(value, self.__ExpectedValue(sensors, error_threshold), 8)

        num_zeroed = sum(1 for sensor in sensors if sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR) == 0.0)
        self.assertGreater(num_zeroed, 0)
        self.assertLess(num_zeroed, len(sensors))

    def test_ErrorThresholdOmittedUsesDefault(self):
        # "error_threshold" is not provided in the sensor parameters at all -> sensor creation
        # should still succeed and fall back to Sensor.DefaultErrorThreshold
        response, sensors = self.__CreateResponseAndSensors(None, "OmittedDefault")

        for sensor in sensors:
            self.assertEqual(sensor.GetErrorThreshold(), KratosSI.Sensors.Sensor.DefaultErrorThreshold)

        value = response.CalculateValue(self.model_part)
        self.assertAlmostEqual(value, self.__ExpectedValue(sensors, KratosSI.Sensors.Sensor.DefaultErrorThreshold), 8)

        # the default threshold is negligibly small, so it should behave exactly like the
        # "below all errors" case: none of the SENSOR_ERROR values get zeroed out.
        for sensor, raw_error in zip(sensors, self.raw_errors):
            self.assertAlmostEqual(sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR), raw_error, 8)

class TestDamageDetectionResponse(kratos_unittest.TestCase):
    def test_DamageResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files/optimization_parameters_p_norm.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.0009799118589232621, 6)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].data

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += delta
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] -= delta

class TestDamageDetectionResponseSensorNormalized(kratos_unittest.TestCase):
    def test_DamageResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files/optimization_parameters_p_norm_local_normalization.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.562499999, 6)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].data

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += delta
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] -= delta

class TestDamageDetectionResponseStrainSensor(kratos_unittest.TestCase):
    def test_DamageResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files_2/optimization_parameters_p_norm.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 2.7829118764552163e-08, 10)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].data

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += delta
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] -= delta

class TestDamageDetectionResponseStrainSensorNormalized(kratos_unittest.TestCase):
    def test_DamageResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files_2/optimization_parameters_p_norm_maximum_normalization.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.140624999, 6)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].data

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += delta
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] -= delta

if __name__ == "__main__":
    kratos_unittest.main()