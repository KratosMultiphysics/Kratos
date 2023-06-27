from abc import ABC, abstractmethod

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
    def Solve(self) -> bool:
        pass

    def GetProcessesOrder(self) -> 'list[str]':
        return ["auxiliary_processes", "output_processes"]

    def CallOnAllProcesses(self, process_types: 'list[str]', *args, **kwargs):
        for process_type in process_types:
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), *args, **kwargs)

