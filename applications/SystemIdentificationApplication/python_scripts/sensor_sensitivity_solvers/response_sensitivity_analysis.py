import abc

import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def CalculateGradient(self, response_function: Kratos.AdjointResponseFunction) -> 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]':
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