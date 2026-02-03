from typing import Optional
import csv
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.average_damage_detection_response import AverageDamageDetectionResponse

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"AverageMembraneCableConstitutiveParametersDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"AverageMembraneCableConstitutiveParametersDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return AverageMembraneCableConstitutiveParametersDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class AverageMembraneCableConstitutiveParametersDetectionResponse(AverageDamageDetectionResponse):
    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS_X, Kratos.YOUNG_MODULUS_Y, Kratos.POISSON_RATIO_XY, Kratos.SHEAR_MODULUS_XY, KratosSM.PRE_STRESS, Kratos.YOUNG_MODULUS, Kratos.POISSON_RATIO, KratosSM.TRUSS_PRESTRESS_PK2]

    def _GetResponsePrefix(self) -> str:
        return "AverageMembraneCableConstitutiveParametersDetectionResponse"