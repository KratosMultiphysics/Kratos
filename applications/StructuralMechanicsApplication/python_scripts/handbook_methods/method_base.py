from abc import ABC, abstractmethod

class HandbookMethod(ABC):
    name: str
    category: str

    @abstractmethod
    def IsApplicable(self, component) -> bool:
        pass

    @abstractmethod
    def Evaluate(self, component):
        pass
