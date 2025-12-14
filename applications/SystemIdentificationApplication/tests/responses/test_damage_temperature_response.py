import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors

class TestDamageTemperatureDetectionAdjointResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part_1 = cls.model.CreateModelPart("sensors1")
        cls.sensor_model_part_2 = cls.model.CreateModelPart("sensors2")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

        parameters1 = [
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
            }""")
        ]

        parameters2 = [
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "weight"       : 2.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "weight"       : 3.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors1 = CreateSensors(cls.sensor_model_part_1, cls.model_part, parameters1)
        cls.sensors2 = CreateSensors(cls.sensor_model_part_2, cls.model_part, parameters2)

        cls.adjoint_response_function1 = KratosSI.Responses.MeasurementResidualResponseFunction(3.0)

        for i, sensor in enumerate(cls.sensors1):
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, i * 15 - 10)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
            cls.adjoint_response_function1.AddSensor(sensor)
        cls.adjoint_response_function1.Initialize()

        cls.adjoint_response_function2 = KratosSI.Responses.MeasurementResidualResponseFunction(4.0)
        for i, sensor in enumerate(cls.sensors2):
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, i * 100 + 20)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
            cls.adjoint_response_function2.AddSensor(sensor)
        cls.adjoint_response_function2.Initialize()

    def setUp(self) -> None:
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
            node.SetSolutionStepValue(Kratos.TEMPERATURE, node.Id * 10.0)
        self.ref_value1 = self.adjoint_response_function1.CalculateValue(self.model_part)
        self.ref_value2 = self.adjoint_response_function2.CalculateValue(self.model_part)

    def test_CalculateValue(self):
        value1 = 0.0
        for sensor in self.sensors1:
            value1 += (0.5 * sensor.GetWeight() * (sensor.CalculateValue(self.model_part) - sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)) ** 2) ** 3.0
        self.assertAlmostEqual(self.ref_value1, value1 ** (1 / 3), 5)

        value2 = 0.0
        for sensor in self.sensors2:
            value2 += (0.5 * sensor.GetWeight() * (sensor.CalculateValue(self.model_part) - sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)) ** 2) ** 4.0
        self.assertAlmostEqual(self.ref_value2, value2 ** (1 / 4), 5)

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
            perturbed_value = self.adjoint_response_function1.CalculateValue(self.model_part)
            global_fd_x_sensitivities[node.Id] = (perturbed_value - self.ref_value1) / delta
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
            perturbed_value = self.adjoint_response_function1.CalculateValue(self.model_part)
            global_fd_y_sensitivities[node.Id] = (perturbed_value - self.ref_value1) / delta
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)

            analytical_x_sensitivities[node.Id] = 0.0
            analytical_y_sensitivities[node.Id] = 0.0

        for element in self.model_part.Elements:
            self.adjoint_response_function1.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                analytical_x_sensitivities[node.Id] += response_sensitivities[i * 6]
                analytical_y_sensitivities[node.Id] += response_sensitivities[i * 6 + 1]

        for node_id in global_fd_x_sensitivities.keys():
            self.assertAlmostEqual(global_fd_x_sensitivities[node_id], analytical_x_sensitivities[node_id], 5)
            self.assertAlmostEqual(global_fd_y_sensitivities[node_id], analytical_y_sensitivities[node_id], 5)

@kratos_unittest.skipIfApplicationsNotAvailable("ConstitutiveLawsApplication")
class TestDamageTemperatureDetectionResponse(kratos_unittest.TestCase):
    def test_DamageTemperatureResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files_7/optimization_parameters_p_norm.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("sum", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 168809.49890615378, 7)

            response.CalculateGradient(var)

            gradients1 = var[Kratos.YOUNG_MODULUS].Evaluate()
            gradients2 = var[Kratos.TEMPERATURE].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = (1 + delta) * E_orig
                sensitivity = ((response.CalculateValue() - ref_value) / ( delta * E_orig))
                self.assertAlmostEqual(gradients1[index], sensitivity, 7)
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

            delta = 1e-3
            for index, node in enumerate(model_part.Nodes):
                node: Kratos.Node
                nodal_temp = node.GetSolutionStepValue(Kratos.TEMPERATURE)
                node.SetSolutionStepValue(Kratos.TEMPERATURE,   nodal_temp * (1 + delta))
                sensitivity = ((response.CalculateValue() - ref_value) / (delta * nodal_temp))
                self.assertAlmostEqual(gradients2[index], sensitivity, 1)
                node.SetSolutionStepValue(Kratos.TEMPERATURE,  nodal_temp)
        

@kratos_unittest.skipIfApplicationsNotAvailable("ConstitutiveLawsApplication")
class TestDamageTemperatureDetectionResponseStrainSensor(kratos_unittest.TestCase):
    def test_DamageTemperatureResponse(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files_8/optimization_parameters_p_norm.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)

            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("sum", ResponseRoutine)

            var = objective.GetRequiredPhysicalGradients()
            response = objective.GetReponse()
            model_part = response.GetInfluencingModelPart()

            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 82.83951393136022, 7)

            response.CalculateGradient(var)

            gradients1 = var[Kratos.YOUNG_MODULUS].Evaluate()
            gradients2 = var[Kratos.TEMPERATURE].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                E_orig = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] = (1 + delta) * E_orig
                sensitivity = ((response.CalculateValue() - ref_value) / (delta * E_orig))
                self.assertAlmostEqual(gradients1[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] = E_orig

            delta = 1e-3
            for index, node in enumerate(model_part.Nodes):
                node: Kratos.Node
                nodal_temp = node.GetSolutionStepValue(Kratos.TEMPERATURE)
                node.SetSolutionStepValue(Kratos.TEMPERATURE,   nodal_temp * (1 + delta))
                sensitivity = ((response.CalculateValue() - ref_value) / (delta * nodal_temp))
                self.assertAlmostEqual(gradients2[index], sensitivity, 3)
                node.SetSolutionStepValue(Kratos.TEMPERATURE,  nodal_temp)

if __name__ == "__main__":
    kratos_unittest.main()