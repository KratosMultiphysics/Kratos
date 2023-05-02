from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunction(ABC):
    """Base response function.

    This reponse function is the base response function. This is assumed to have following responsibilities.
        1. CalculateValue for a fresh design. (@see CalculateValue)
        2. CalculateSensitivity for a fresh design (@see CalculateSensitivity)

    This control should only work on one model part. Hence, if multiple model parts required then,
    a single model part should be created using Kratos.ModelPartOperationUtilities.
    """

    @abstractmethod
    def Initialize(self) -> None:
        """Initializes the response.

        This method initializes the response. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def Check(self) -> None:
        """Checks the response.

        This method checks the response. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def Finalize(self) -> None:
        """Finalizes the response.

        This method finalizes the response. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def GetDependentPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        """Returns all the dependent physical kratos variables of the response.

        This method should return all the dependent physical variables of the response to make sure
        that sensitivities w.r.t. those variables are requested from the response. If sensitivities
        w.r.t. other variables are required from this response, they are ASSUMED TO BE ZERO.

        Please return all the dependent physical variables, eventhough the gradient computation is not yet implemented
        for some to avoid future bugs.

        Returns:
            list[SupportedSensitivityFieldVariableTypes]: All dependent physical variables of the response.
        """
        pass

    @abstractmethod
    def CalculateValue(self) -> float:
        """Calculates the response value.

        This method should always calculate the response value assuming the domain has changed.

        Returns:
            float: Calculated response value.
        """
        pass

    @abstractmethod
    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]') -> None:
        """Calculate gradient w.r.t. given physical variables.

        This method should always calculate the sensitivities w.r.t. requested physical variables on the given container expressions
        in the collective expression. An error should be thrown if sensitivitiy computation w.r.t. one or more given physical variables
        are not implemented.

        This method should always calculate the sensitivities assuming the domain has changed.

        physical_variable_collective_expressions is a map of physical variables, and their domains. The domains are represented by a CollectiveExpression
        which contains list of empty ContainerExpression. Each empty ContainerExpression contains details of the model part's nodes/conditions/element/properties
        container for which the sensitivities w.r.t. physical variable requested.

        Args:
            physical_variable_collective_expressions (dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]): Output containing calculated sensitivities w.r.t. requested physical variables.
        """
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        """Returns the model part for which this response is computed on.

        Returns:
            Kratos.ModelPart: Response function model part.
        """
        pass