from abc import ABC, abstractmethod

import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class Control(ABC):
    def Initialize(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    @abstractmethod
    def GetRequiredSensitivityFieldVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        pass

    @abstractmethod
    def GetCollectiveExpression(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        pass

    @abstractmethod
    def MapFirstDerivative(self, sensitivity_variable_first_derivative_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', mapped_first_derivative: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        pass

    @abstractmethod
    def Update(self, update: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        pass

    def _CheckSensitivityCollectiveExpressions(self, sensitivity_variable_first_derivative_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        required_sensitivity_variables = self.GetRequiredSensitivityFieldVariables()
        required_collective_expression = self.GetCollectiveExpression()

        # check whether correct sensitivity variable sensitivities are provided
        if sorted(self.GetRequiredSensitivityFieldVariables(), key=lambda x:x.Name()) != sorted(sensitivity_variable_first_derivative_map.keys(), key=lambda x:x.Name()):
            raise RuntimeError("Invalid sensitivity variable collective expressions provided. Required sensitivity variables:\n\t" + "\n\t".join(v.Name() for v in required_sensitivity_variables) + "\nProvided sensitivity variables:\n\t" + "\n\t".join([v.Name() for v in sensitivity_variable_first_derivative_map.keys()]))

        for sensitivity_collective_expression in sensitivity_variable_first_derivative_map.values():
            # now check sensitivities for all required model parts are provided in the same
            if len(required_collective_expression.GetContainerExpressions()) != len(sensitivity_collective_expression.GetContainerExpressions()):
                raise RuntimeError(f"Invalid sensitivity collective expressions. Required sensitivity expressions:\n{required_collective_expression}\nProvided sensitivity model parts:\n{sensitivity_collective_expression}")

            for required_container_expression, provided_container_expression in zip(required_collective_expression.GetContainerExpressions(), sensitivity_collective_expression.GetContainerExpressions()):
                if required_container_expression.GetModelPart() != provided_container_expression.GetModelPart():
                    raise RuntimeError(f"Invalid sensitivity collective expressions. Required sensitivity expressions:\n{required_collective_expression}\nProvided sensitivity model parts:\n{sensitivity_collective_expression}")

    def _CheckUpdateCollectiveExpressions(self, update: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        required_collective_expression = self.GetCollectiveExpression()

        # now check update for all required model parts are provided in the same
        # order as in the GetModelParts method
        if len(required_collective_expression.GetContainerExpressions()) != len(update.GetContainerExpressions()):
            raise RuntimeError(f"Invalid update collective expressions. Required update expressions:\n{required_collective_expression}\nProvided update model parts:\n{update}")

        for required_container_expression, provided_container_expression in zip(required_collective_expression.GetContainerExpressions(), update.GetContainerExpressions()):
            if required_container_expression.GetModelPart() != provided_container_expression.GetModelPart():
                raise RuntimeError(f"Invalid update collective expressions. Required update expressions:\n{required_collective_expression}\nProvided update model parts:\n{update}")

