# STD Imports
from typing import Optional
from abc import ABC, abstractmethod

# Kratos Imports
import KratosMultiphysics


class StepController(ABC):
    def __init__(self, parameters: KratosMultiphysics.Parameters) -> None:
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

    @abstractmethod
    def GetNextStep(self,
                    begin: float,
                    end: float,
                    converged: bool,
                    substep_index: int) -> Optional[float]:
        pass

    @classmethod
    @abstractmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters()


class DefaultStepController(StepController):
    def __init__(self, parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__(parameters)
        self.__end_time: float = parameters["end_time"].GetDouble()

    def GetNextStep(self,
                    begin: float,
                    end: float,
                    converged: bool,
                    substep_index: int) -> Optional[float]:
        return end if self.__end_time < end else None

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "end_time" : 0
        }""")


class GeometricStepController(StepController):
    def __init__(self, parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__(parameters)
        self.__end_time: float = parameters["end_time"].GetDouble()
        self.__min_step_size: float = parameters["min_step_size"].GetDouble()
        self.__increase_factor: float = parameters["increase_factor"].GetDouble()
        self.__decrease_factor: float = parameters["decrease_factor"].GetDouble()
        self.__max_substeps: int = parameters["max_substeps"].GetInt()
        self.__last_step_size: Optional[float] = None
        self.__Check()

    def GetNextStep(self,
                    begin: float,
                    end: float,
                    converged: bool,
                    substep_index: int) -> Optional[float]:
        if self.__end_time <= begin or self.__max_substeps <= substep_index:
            return None

        max_step_size: float = end - begin

        if max_step_size < 0:
            KratosMultiphysics.Logger.PrintWarning(
                f"Negative time step from {begin} to {end}.",
                label = type(self).__name__)

        if max_step_size < self.__min_step_size:
            KratosMultiphysics.Logger.PrintWarning(
                f"Time step from {begin} to {end} is smaller than the minimum requested step size {self.__min_step_size}.",
                label = type(self).__name__)

        if not isinstance(self.__last_step_size, float):
            self.__last_step_size = max_step_size

        step_size: float = self.__last_step_size
        if converged:
            step_size *= self.__increase_factor
        else:
            step_size = max(self.__min_step_size, self.__decrease_factor * step_size)

        time: float = min(begin + step_size, end)
        self.__last_step_size = time - begin
        return time

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "end_time" : 0,
            "min_step_size" : 0,
            "increase_factor" : 2e0,
            "decrease_factor" : 5e-1,
            "max_substeps" : 1e1
        }""")

    def __Check(self) -> None:
        if self.__min_step_size < 0:
            raise ValueError(f"\"min_step_size\" ({self.__min_step_size}) must be non-negative.")

        if self.__increase_factor <= 1.0:
            raise ValueError(f"\"increase_factor\" must be greater than 1.")

        if self.__decrease_factor <= 0.0 or 1.0 <= self.__decrease_factor:
            raise ValueError(f"\"decrease_factor\" must be between 0 and 1.")

        if self.__max_substeps < 0:
            raise ValueError(f"\"max_substeps\" must be a non-negative integer.")
