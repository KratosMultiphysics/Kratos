import numpy as np

from abc import ABC, abstractmethod
from typing import Union
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DistanceResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DistanceResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return OneMinusDistanceSquaredResponseFunction(parameters["name"].GetString(), model, parameters["settings"])

class OneMinusDistanceSquaredResponseFunction(ResponseFunction):
    """Base response function.

    This response function is the base response function. This is assumed to have following responsibilities.
        1. CalculateValue for a new design. (@see CalculateValue)
        2. CalculateSensitivity for a new design (@see CalculateSensitivity)

    This response should only work on one model part. Hence, if multiple model parts required then,
    a single model part should be created using Kratos.ModelPartOperationUtilities.
    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{}""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.model_part = self.model["Structure"]

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def CalculateValue(self) -> float:
        """Calculates the response value.

        This method should always calculate the response value assuming the domain has changed.

        Returns:
            float: Calculated response value.
        """
        x0 = self.model_part.GetNode(1).X
        y0 = self.model_part.GetNode(1).Y
        # z0 = self.model_part.GetNode(1).Z
        x1 = self.model_part.GetNode(2).X
        y1 = self.model_part.GetNode(2).Y
        # z1 = self.model_part.GetNode(2).Z
        distance = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        f = (1-distance)**2
        return float(f)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        """Calculate gradient w.r.t. given physical variables.

        This method should always calculate the sensitivities w.r.t. requested physical variables on the given container expressions
        in the collective expression. An error should be thrown if sensitivitiy computation w.r.t. one or more given physical variables
        are not implemented.

        This method should always calculate the sensitivities assuming the domain has changed.

        physical_variable_collective_expressions is a map of physical variables, and their domains. The domains are represented by a CollectiveExpression
        which contains list of empty ContainerExpression. Each empty ContainerExpression contains details of the model part's nodes/conditions/element/properties
        container for which the sensitivities w.r.t. physical variable requested.

        Args:
            physical_variable_collective_expressions (dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]): Output containing calculated sensitivities w.r.t. requested physical variables.
        """
        x0 = self.model_part.GetNode(1).X
        y0 = self.model_part.GetNode(1).Y
        # z0 = self.model_part.GetNode(1).Z
        x1 = self.model_part.GetNode(2).X
        y1 = self.model_part.GetNode(2).Y
        # z1 = self.model_part.GetNode(2).Z
        distance = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        dfdx1 = -2*(1-distance)/distance * x1
        dfdy1 = -2*(1-distance)/distance * y1
        # dfdz1 = -2*(1-distance)/distance * z1
        grad_list = np.array([[0.0, 0.0, 0.0], [ dfdx1, dfdy1, 0.0]], dtype=np.float64)
        expression = physical_variable_collective_expressions[KratosOA.SHAPE].GetContainerExpressions()[0]
        Kratos.Expression.CArrayExpressionIO.Read(expression, grad_list)

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        """Returns the model part which influences the computation of the response value.

        Following two cases are considered:
            1. Responses without adjoint system solve: The value of the response can only be influenced by
               changing the quantities in the evaluated model part (evaluated model part is the one which the
               response value is computed.) Therefore, in this case, this method should return the
               evaluated model part.
            2. Responses with adjoint system solve: The value of the response can be influenced by
                changing quantities in the adjoint/primal model part (Evaluated model part needs to have
                intersection with the adjoint/primal model part). Therefore, in this case, this method should return the
                adjoint analysis model part.

        Returns:
            Kratos.ModelPart: Response function model part which influences the response value.
        """
        return self.model_part


