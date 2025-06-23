from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class Algorithm(ABC):
    """Base class for algorithm.

    Algorithm works purely in the control space. Hence, the ResponseRoutine which
    the algorithm is associated with should convert all physical space gradients to
    control space and control space designs to physical space designs.

    The purpose of the algorithm is to solve a minimization problem. Hence, the used
    ResponseRoutine should standardize the values and gradients it computes accordingly.

    Algorithm should create one MasterControl from all the require controls, and share
    it among all the ResponseRoutines.

    Once it received the objective/constraint values and their gradients, the algorithm
    will compute the new design and request new objective/constraints values and their
    gradients until the specified convergence criteria is met.
    """
    def __init__(self, optimization_problem: OptimizationProblem) -> None:
        self._optimization_problem = optimization_problem

    @abstractmethod
    def Initialize(self) -> None:
        """Initializes the algorithm, ResponseRoutines and MasterControl.
        """
        pass

    @abstractmethod
    def Check(self) -> None:
        """Checks algorithm, ResponseRoutines and MasterControl.
        """
        pass

    @abstractmethod
    def Finalize(self) -> None:
        """Finalizes algorithm ResponseRoutines and MasterControl.
        """
        pass

    @abstractmethod
    def Solve(self) -> bool:
        """Solves the optimization problem.

        This method will have a loop which iterates until the specified
        convergence for the given optimization problem.

        Returns:
            bool: Convergence status. True for converged, False for non-converged.
        """
        pass

    def GetProcessesOrder(self) -> 'list[str]':
        """The order of execution of the process categories.

        This defines the order of execution of the processes defined under "processes"
        in either "kratos_process" or "optimization_data_processes".

        Returns:
            list[str]: List of strings for process categories.
        """
        return ["auxiliary_processes", "output_processes"]

    def CallOnAllProcesses(self, process_types: 'list[str]', *args, **kwargs):
        for process_type in process_types:
            CallOnAll(self._optimization_problem.GetListOfProcesses(process_type), *args, **kwargs)

    def _InitializeIteration(self) -> None:
        for process_type in self._optimization_problem.GetAvailableProcessTypes():
            list(map(lambda x: x.ExecuteInitializeSolutionStep(), self._optimization_problem.GetListOfProcesses(process_type)))

    def _FinalizeIteration(self) -> None:
        for process_type in self._optimization_problem.GetAvailableProcessTypes():
            list(map(lambda x: x.ExecuteFinalizeSolutionStep(), self._optimization_problem.GetListOfProcesses(process_type)))
