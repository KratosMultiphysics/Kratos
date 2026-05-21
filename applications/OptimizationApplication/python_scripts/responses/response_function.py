from abc import ABC, abstractmethod
from typing import Union
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunction(ABC):
    """Base response function.

    This reponse function is the base response function. This is assumed to have following responsibilities.
        1. CalculateValue for a new design. (@see CalculateValue)
        2. CalculateSensitivity for a new design (@see CalculateSensitivity)

    This response should only work on one model part. Hence, if multiple model parts required then,
    a single model part should be created using Kratos.ModelPartOperationUtilities.
    """
    def __init__(self, response_name: str) -> None:
        self.__name = response_name

    def GetName(self) -> str:
        return self.__name

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
    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
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
        pass

    @abstractmethod
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
        pass

    def GetChildResponses(self) -> 'list[ResponseFunction]':
        """Returns the list of child responses.

        This method returns list of child responses. If the current response is a leaf response,
        then it will return empty list.

        Needs to be implemented in non-leaf response functions.

        Returns:
            list[ResponseFunction]: List of child responses.
        """
        return []
