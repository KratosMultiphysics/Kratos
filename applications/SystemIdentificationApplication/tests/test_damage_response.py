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
            with open("auxiliary_files/system_identification/optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            
            analysis.Initialize()
            analysis.Check()
            objective: ResponseRoutine = analysis.optimization_problem.GetComponent("damage_response", ResponseRoutine)
            var = objective.GetRequiredPhysicalGradients()
            print(var)
            response = analysis.optimization_problem.GetResponse("damage_response")
            ref_value = response.CalculateValue()
            self.assertAlmostEqual(ref_value, 0.0007924977797682586, 12)
            sensitivity = analysis.optimization_problem.GetComponent("master_control", MasterControl).GetEmptyField()
            response.CalculateGradient(var)
            gradients = var[Kratos.YOUNG_MODULUS].Evaluate()

            model_part = response.GetInfluencingModelPart()

            for index, element in enumerate(model_part.Elements):
                print(element.Properties[Kratos.YOUNG_MODULUS])
                element.Properties[Kratos.YOUNG_MODULUS] += 1e-6
                print(element.Properties[Kratos.YOUNG_MODULUS])
                new_value = response.CalculateValue()
                print(new_value, ref_value)
                sensitivity = (new_value - ref_value) / 1e-6
                print(sensitivity)
                print(gradients[index])
                raise RuntimeError(1)
                element.Properties[Kratos.YOUNG_MODULUS] -= 1
                print(element.Properties[Kratos.YOUNG_MODULUS])
                

            # for 


            
            
    


if __name__ == "__main__":
    kratos_unittest.main()