import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

import time as timer
import numpy as np
from typing import Optional

def Factory(model: KM.Model, parameters: KM.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"MassResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MassResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return DummyResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class DummyResponseFunction(ResponseFunction):

    def __init__(self, name: str, model: KM.Model, parameters: KM.Parameters):

        super().__init__(name)

        default_settings = KM.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for MassResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[KM.ModelPart] = None        

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KOA.SHAPE]

    def GetEvaluatedModelPart(self) -> KM.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call MassResponseFunction::Initialize first.")
        return self.model_part

    def Initialize(self):
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        if self.model_part is None:
            raise RuntimeError("Please call MassResponseFunction::Initialize first.")

    def Finalize(self) -> None:
        pass

    def CalculateValue(self):
        Logger.PrintInfo("MassResponseFunction:CalculateValue: Starting value calculation for response ", self.GetName())
        startTime = timer.time()
        Logger.PrintInfo("MassResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")        
        return 0.0
    
    def GetAnalysisModelPart(self) -> None:
        return None

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KOA.CollectiveExpression]') -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, False)

        # calculate the gradients
        for physical_variable, merged_model_part in merged_model_part_map.items():
            expression_list = physical_variable_collective_expressions[physical_variable].GetContainerExpressions()
            
            for expression in expression_list:
                num_nodes = expression.GetModelPart().NumberOfNodes()
                gradient_array = np.zeros((num_nodes, 3))
                KM.Expression.CArrayExpressionIO.Read(expression, gradient_array)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"