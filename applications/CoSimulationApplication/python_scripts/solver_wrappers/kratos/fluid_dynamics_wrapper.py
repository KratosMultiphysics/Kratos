# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    raise ImportError("The FluidDynamicsApplication is not available!")
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities

def Create(settings, model, solver_name):
    return FluidDynamicsWrapper(settings, model, solver_name)

class FluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the FluidDynamicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        if not KM.IsDistributedRun():
            return KM.ParallelEnvironment.GetDataCommunicator("Serial")

        # now we know that Kratos runs in MPI
        parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        solver_wrapper_settings = self.settings["solver_wrapper_settings"]

        # first check if the solver uses MPI
        if parallel_type != "MPI":
            if solver_wrapper_settings.Has("data_communicator_creation"):
                raise Exception('"data_communicator_creation" was specified but solver is not running in MPI!')
            return data_communicator_utilities.GetRankZeroDataCommunicator()

        # now we know that the solver uses MPI, only question left is whether to use all ranks or a subset
        if self.project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid":
            model_import_settings = self.project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]
        else:
            model_import_settings = self.project_parameters["solver_settings"]["model_import_settings"]

        self._CheckDataCommunicatorIsConsistentlyDefined(model_import_settings, solver_wrapper_settings)

        if solver_wrapper_settings.Has("data_communicator_creation"):
            return data_communicator_utilities.CreateDataCommunicatorWithNProcesses(solver_wrapper_settings["data_communicator_creation"])
        else:
            return KM.ParallelEnvironment.GetDataCommunicator("World")
