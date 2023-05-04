from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class Algorithm(ABC):
    def __init__(self, optimization_problem: OptimizationProblem) -> None:
        self._optimization_problem = optimization_problem

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
    def SolveOptimizationProblem(self) -> None:
        pass

    def GetProcessesOrder(self) -> 'list[str]':
        return ["auxiliary_processes", "output_processes"]

    def ExecuteBeforeSolutionLoop(self):
        for process_type in self.GetProcessesOrder():
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeSolutionLoop)

    def ExecuteInitializeSolutionStep(self):
        for process_type in self.GetProcessesOrder():
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitializeSolutionStep)

    def ExecuteBeforeOutputStep(self):
        for process_type in self.GetProcessesOrder():
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeOutputStep)

    def ExecuteBeforeOutputStep(self):
        for process_type in self.GetProcessesOrder():
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeOutputStep)