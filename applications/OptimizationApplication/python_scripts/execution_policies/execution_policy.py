from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class ExecutionPolicy(ABC):
    """Base class for execution policy

    This represents the base class for the execution policy which
    can be used to execute primal analysis to obtain response function
    values and their gradients.

    """
    def __init__(self, execution_policy_name: str) -> None:
        self.__name = execution_policy_name

    def GetName(self) -> str:
        return self.__name

    @abstractmethod
    def Initialize(self) -> None:
        pass

    @abstractmethod
    def Check(self) -> None:
        pass

    @abstractmethod
    def Finalize(self) -> None:
        pass

    @abstractmethod
    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        """Returns the analysis model part.

        Each execution policy will be using a domain to compute the primal model part.
        This method should return that domain as a Kratos::ModelPart.

        Returns:
            Kratos.ModelPart: Domain used to compute the primal solution.
        """
        pass

    @abstractmethod
    def Execute(self) -> None:
        """Solves the associated primal problem.
        """
        pass