import abc

import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def CalculateGradient(self, response_function: Kratos.AdjointResponseFunction) -> None:
        """Calculate the gradient using the provided adjoint response function.

        This method calculate all the sensitivities using the given adjoint response function.

        Args:
            response_function (Kratos.AdjointResponseFunction): Adjoint response function to be used for sensitivity computation
        """
        pass

    def GetGradient(self, sensitivity_variable: SupportedSensitivityFieldVariableTypes, gradient_expression: ContainerExpressionTypes) -> None:
        """Returns the gradients in the domain represented by gradient_expression container expression.

        Args:
            sensitivity_variable (SupportedSensitivityFieldVariableTypes): Sensitivity variable
            gradient_expression (ContainerExpressionTypes): Container expression to hold the gradients.
        """
        if isinstance(gradient_expression, Kratos.Expression.NodalExpression):
            Kratos.Expression.VariableExpressionIO.Read(gradient_expression, sensitivity_variable, True)
        else:
            Kratos.Expression.VariableExpressionIO.Read(gradient_expression, sensitivity_variable)

    def GetProcessesOrder(self) -> 'list[str]':
        """The order of execution of the process categories.

        This defines the order of execution of the processes defined under "processes"
        in either "kratos_process" or "optimization_data_processes".

        Returns:
            list[str]: List of strings for process categories.
        """
        return ["auxiliary_processes", "output_processes"]