from abc import ABC
from abc import abstractmethod

class Algorithm(ABC):
    def __init__(self):
        self.__list_of_constraints = []
        self.__list_of_objectives = []
        self.__list_of_controllers = []
        self.__name = ""

    def AddVariables(self):
        pass

    def AddDofs(self):
        pass

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def SetName(self, name: str):
        self.__name = name

    def SetObjectives(self, list_of_objectives):
        self.__list_of_objectives = list_of_objectives

    def SetConstraints(self, list_of_constraints):
        self.__list_of_constraints = list_of_constraints

    def SetControllers(self, list_of_controllers):
        self.__list_of_controllers = list_of_controllers

    def GetName(self) -> str:
        return self.__name

    def GetObjectives(self):
        return self.__list_of_objectives

    def GetConstraints(self):
        return self.__list_of_constraints

    def GetControllers(self):
        return self.__list_of_controllers

    @abstractmethod
    def GetMinimumBufferSize(self) -> int:
        pass

    @abstractmethod
    def Check(self):
        pass

    @abstractmethod
    def SolveSolutionStep(self):
        pass

    @abstractmethod
    def IsConverged(self) -> bool:
        pass