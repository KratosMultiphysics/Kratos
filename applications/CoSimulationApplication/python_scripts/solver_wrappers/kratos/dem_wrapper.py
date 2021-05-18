# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing StructuralMechanics
if not CheckIfApplicationsAvailable("DEMApplication"):
    raise ImportError("The DEMApplication is not available!")

from KratosMultiphysics import DEMApplication
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

def Create(settings, model, solver_name):
    return DEMWrapper(settings, model, solver_name)

class DEMWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the DEMApplication of Kratos"""

    def _CreateAnalysisStage(self):
        dem_analysis_module = DEMAnalysisStage

        if self.settings["solver_wrapper_settings"].Has("working_directory"):
            working_dir = self.settings["solver_wrapper_settings"]["working_directory"].GetString()

            class DEMAnalysisStageWithWorkingDir(DEMAnalysisStage):
                @classmethod
                def GetMainPath(self):
                    return working_dir

            dem_analysis_module = DEMAnalysisStageWithWorkingDir

        return dem_analysis_module(self.model, self.project_parameters)

    def Initialize(self):
        super().Initialize()

        # save nodes in model parts which need to be moved while simulating
        self.list_of_nodes_in_move_mesh_model_parts = [self.model[mp_name].Nodes for mp_name in self.settings["solver_wrapper_settings"]["move_mesh_model_part"].GetStringArray()]


    def SolveSolutionStep(self):
        # move the rigid wall object in the dem mp w.r.t. the current displacement and velocities
        for model_part_nodes in self.list_of_nodes_in_move_mesh_model_parts:
            DEMApplication.MoveMeshUtility().MoveDemMesh(model_part_nodes,True)
        # solve DEM
        super().SolveSolutionStep()

