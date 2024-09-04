from abc import ABC, abstractmethod
import numpy
from math import sqrt

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors


class TestDamageDetectionAdjointResponseFunction(kratos_unittest.TestCase):
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

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = GetSensors(cls.model_part, parameters)

        cls.adjoint_response_function = KratosSI.Sensors.MeasurementResidualResponseFunction(3.0)

        for i, sensor in enumerate(cls.sensors):
            sensor.SetValue(KratosSI.SENSOR_MEASURED_VALUE, i * 15 - 10)
            cls.adjoint_response_function.AddSensor(sensor)

        cls.adjoint_response_function.Initialize()

    def setUp(self) -> None:
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
        self.ref_value = self.adjoint_response_function.CalculateValue(self.model_part)

    def test_CalculateValue(self):
        value = 0.0
        for sensor in self.sensors:
            value += (0.5 * (sensor.CalculateValue(self.model_part) - sensor.GetValue(KratosSI.SENSOR_MEASURED_VALUE)) ** 2) ** 3.0
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
            self.assertAlmostEqual(ref_value, 0.017778076997606067, 6)

            response.CalculateGradient(var)

            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            delta = 1e-8
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += delta
                sensitivity = ((response.CalculateValue() - ref_value) / delta)
                self.assertAlmostEqual(gradients[index], sensitivity, 6)
                element.Properties[Kratos.YOUNG_MODULUS] -= delta

if __name__ == "__main__":
    kratos_unittest.main()