from abc import ABC, abstractmethod
import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl



class TestDamageDetectionResponseBase(kratos_unittest.TestCase, ABC):
    def test_damage_response(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files/optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            
            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)
            var = objective.GetRequiredPhysicalGradients()
            response = analysis.optimization_problem.GetResponse("damage_response")
            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.017887292422409877, 6)
            sensitivity = analysis.optimization_problem.GetComponent("master_control", MasterControl).GetEmptyField()
            response.CalculateGradient(var)
            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()
            model_part = response.GetInfluencingModelPart()
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] = 20.0 + 1e-6
                new_value = response.CalculateValue()
                print(new_value, ref_value)
                sensitivity = (new_value - ref_value) / 1e-6
                print(sensitivity, gradients[index])
                # self.assertAlmostEqual(gradients[index], sensitivity, 12)
                element.Properties[Kratos.YOUNG_MODULUS] = 20.0

    def test_damage_response_p_norm(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files/optimization_parameters_p_norm.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            
            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)
            var = objective.GetRequiredPhysicalGradients()
            response = analysis.optimization_problem.GetResponse("damage_response")
            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.017887292422409877, 12)
            sensitivity = analysis.optimization_problem.GetComponent("master_control", MasterControl).GetEmptyField()
            response.CalculateGradient(var)
            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()
            model_part = response.GetInfluencingModelPart()
            for index, element in enumerate(model_part.Elements):
                element.Properties[Kratos.YOUNG_MODULUS] += 1e-6
                new_value = response.CalculateValue()
                sensitivity = (new_value - ref_value) / 1e-6
                # self.assertAlmostEqual(gradients[index], sensitivity, 12)
                element.Properties[Kratos.YOUNG_MODULUS] -= 1e-6

if __name__ == "__main__":
    kratos_unittest.main()