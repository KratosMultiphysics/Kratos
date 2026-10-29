from abc import ABC, abstractmethod

class VerificationMethod(ABC):
    name: str
    category: str

    def CheckCompatibility(self, component) -> None:
        """Raise if the method cannot be used for this component (configuration error).
        Called once before the solution loop.
        """

    @abstractmethod
    def IsApplicable(self, component) -> bool:
        """False if the method does not apply to the current load state (result: RF = inf)."""
        pass

    @abstractmethod
    def Evaluate(self, component):
        pass
