import abc

import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def CalculateGradient(self, response_function: Kratos.AdjointResponseFunction) -> None:
        """Returns the gradient computed by using the provided adjoint response function.

        This method returns all the sensitivities computed from the given adjoint response function. The return
        will be a dictionary having a pair with the variable as the key and a ContainerExpression as the value.
        ContainerExpression will carry the sensitivity values.

        Args:
            response_function (Kratos.AdjointResponseFunction): Adjoint response function to be used for sensitivity computation

        Returns:
            dict[str, ContainerExpressionTypes]: Sensitivities dictionary with variable and sensitivities as the pair.
        """
        pass

    @abc.abstractmethod
    def GetGradient(self, sensitivity_variable: SupportedSensitivityFieldVariableTypes, container_expression: ContainerExpressionTypes) -> None:
        pass

    def GetProcessesOrder(self) -> 'list[str]':
        """The order of execution of the process categories.

        This defines the order of execution of the processes defined under "processes"
        in either "kratos_process" or "optimization_data_processes".

        Returns:
            list[str]: List of strings for process categories.
        """
        return ["auxiliary_processes", "output_processes"]