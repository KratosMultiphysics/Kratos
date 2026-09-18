from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

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
    # Whether this algorithm's Solve() calls SaveRestartFile(). False by default; the three
    # gradient-descent-family algorithms whose Solve() does call it set it True, and it is
    # inherited unchanged by every subclass that only overrides ComputeControlUpdate (Adam,
    # Nesterov, momentum-relaxed-gradient-projection). OptimizationAnalysis._CreateRestart() checks
    # this to fail loudly rather than silently never writing a checkpoint.
    SUPPORTS_RESTART = False

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

    def SaveRestartFile(self) -> None:
        """Writes a restart checkpoint if save_restart is on and this step is a save step.

        Must be called right after UpdateControl(), never before: Output() writes "control_field"
        before that step's delta is applied, so the checkpoint is only self-consistent (mesh and
        BufferedDict agreeing on the same design, no pending delta) once UpdateControl() has run.
        """
        restart_output = self._optimization_problem.GetRestartOutput()
        if restart_output is None or not restart_output.IsSaveStep():
            return

        # UpdateControl() calls StoreData(perform_store_data_recursively=False), which only pushes
        # the just-applied delta into each control's own TensorAdaptor, not all the way down to the
        # mesh -- normal operation never notices, since every later read goes through that same
        # TensorAdaptor, never the mesh directly. This checkpoint reads the mesh directly (it saves
        # the whole Model), so push through Control.Update() first -- it reads the parts
        # (GetTensorAdaptors()), already correct after UpdateControl(), not .data, and each
        # control's own Update() writes them to the mesh.
        self.master_control.Update(self.GetCurrentControlField())

        # refresh "control_field" to the just-applied design; Output() wrote the pre-update one.
        ComponentDataView("algorithm", self._optimization_problem).GetBufferedData().SetValue(
            "control_field",
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self.GetCurrentControlField()),
            overwrite=True)

        restart_output.SaveRestartFile()
