import abc
import typing

class ConvergenceCriterion(abc.ABC):
    """Initializes the convergence criteria
    """
    @abc.abstractmethod
    def Initialize(self) -> None:
        pass

    """This method returns whether the convergence criteria detected the optimization
    has converged.
    """
    @abc.abstractmethod
    def IsConverged(self) -> bool:
        pass

    """Finalizes the convergence criteria """
    @abc.abstractmethod
    def Finalize(self) -> None:
        pass

    """Returns a list[(str, Any)] with information about the current convergence."""
    @abc.abstractmethod
    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        pass