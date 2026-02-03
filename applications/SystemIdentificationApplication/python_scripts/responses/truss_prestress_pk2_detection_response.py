import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"TrussPrestressPK2DetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"TrussPrestressPK2DetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return TrussPrestressPK2DetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class TrussPrestressPK2DetectionResponse(DamageDetectionResponse):
    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSM.TRUSS_PRESTRESS_PK2]

    def _GetResponsePrefix(self) -> str:
        return "TrussPrestressPK2DetectionResponse"