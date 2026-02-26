import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DamageTemperatureDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DamageTemperatureDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return DamageTemperatureDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class DamageTemperatureDetectionResponse(DamageDetectionResponse):
    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS, Kratos.TEMPERATURE]

    def _GetResponsePrefix(self) -> str:
        return "DamageTemperatureDetectionResponse"