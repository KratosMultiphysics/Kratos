# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing StructuralMechanics
if not CheckIfApplicationsAvailable("GeoMechanicsApplication"):
    raise ImportError("The GeoMechanicsApplication is not available!")
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities

def Create(settings, model, solver_name):
    return GeoMechanicsWrapper(settings, model, solver_name)

class GeoMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the GeoMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return GeoMechanicsAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        if not KM.IsDistributedRun():
            return KM.ParallelEnvironment.GetDataCommunicator("Serial")

        # check if the solver uses MPI
        parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        if parallel_type == "MPI":
            raise RuntimeError("MPI is not supported (yet) by the GeoMechanicsApplication.")

        return super()._GetDataCommunicator()
