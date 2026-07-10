import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"TemperatureDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"TemperatureDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return TemperatureDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class TemperatureDetectionResponse(DamageDetectionResponse):
    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.TEMPERATURE]

    def _GetResponsePrefix(self) -> str:
        return "TemperatureDetectionResponse"