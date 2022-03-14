# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing StructuralMechanics
if not CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    raise ImportError("The StructuralMechanicsApplication is not available!")
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities

def Create(settings, model, solver_name):
    return StructuralMechanicsWrapper(settings, model, solver_name)

class StructuralMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the StructuralMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return StructuralMechanicsAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        if not KM.IsDistributedRun():
            return KM.ParallelEnvironment.GetDataCommunicator("Serial")

        # now we know that Kratos runs in MPI
        parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # first check if the solver uses MPI
        if parallel_type != "MPI":
            return data_communicator_utilities.GetRankZeroDataCommunicator()

        # now we know that the solver uses MPI, only question left is whether to use all ranks or a subset
        self._CheckDataCommunicatorIsConsistentlyDefined(self.project_parameters["solver_settings"]["model_import_settings"], self.settings["mpi_settings"])

        return super()._GetDataCommunicator()
