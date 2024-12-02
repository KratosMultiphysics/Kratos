"""
This response is a duplicate of the face_angle response in ShapeOptApplication.
Created on 20.11.2024
"""
from typing import Optional
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"TargetGeometryResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"TargetGeometryResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return TargetGeometryResponse(parameters["name"].GetString(), model, parameters["settings"])

class TargetGeometryResponse(ResponseFunction):
    """
    Optimize towrds a target geometry.

    """

    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_name"       : "UNKNOWN_NAME"
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.response_settings = parameters
        self.model = model

        self.value = None
        self.target_geometry = None
        self._model_part_name = parameters["model_part_name"].GetString()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self._model_part_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]
    
    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call TargetGeometryResponse::Initialize first.")
        return self.model_part
    
    def GetAnalysisModelPart(self) -> None:
        return None
    
    def CreateTargetGeometry(self) -> None:
        self.target_geometry = np.zeros((self.model_part.NumberOfNodes(),3), dtype=float)
        nodes = self.model_part.GetNodes()
        for i, node in enumerate(nodes):
            self.target_geometry[i,0] = node.X * 0.75
            self.target_geometry[i,1] = node.Y * 2
            self.target_geometry[i,2] = node.Z
    
    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        
        self.CreateTargetGeometry()

        # self.response_function_utility = KratosOA.ResponseUtils.TargetGeometryResponseUtils(only_part, self.response_settings)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def CalculateValue(self) -> float:
        self.value = 0.0
        nodes = self.model_part.GetNodes()
        target_geometry = self.target_geometry
        for i, node in enumerate(nodes):
            node_coords = np.array([node.X, node.Y, node.Z])
            target_coords = target_geometry[i,:]
            self.value += np.linalg.norm(target_coords - node_coords)**2

        print("CalculateValue :: response value ", self.value)
        return self.value
    
    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        print("target_geometry_response.py :: CalculateGradient")

        target_geometry = self.target_geometry
        gradient_field = np.zeros(target_geometry.shape, dtype=float)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.SHAPE_SENSITIVITY, self.model_part.Nodes)

        nodes = self.model_part.GetNodes()
        for i, node in enumerate(nodes):
            node_coords = np.array([node.X, node.Y, node.Z])
            gradient_field[i,:] = - 2 * (target_geometry[i,:] - node_coords)
        
        # import pdb
        # pdb.set_trace()
        
        # Kratos.VariableUtils().SetNonHistoricalVariable(Kratos.SHAPE_SENSITIVITY, gradient_field, self.model_part.Nodes)
        
        physical_variable = KratosOA.SHAPE
        for container_expression in physical_variable_collective_expressions[physical_variable].GetContainerExpressions():
            Kratos.Expression.CArrayExpressionIO.Read(container_expression, gradient_field)


                