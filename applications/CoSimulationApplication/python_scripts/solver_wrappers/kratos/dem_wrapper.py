from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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

def Create(settings, solver_name):
    return DEMWrapper(settings, solver_name)

class DEMWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the DEMApplication of Kratos"""

    def _CreateAnalysisStage(self):
        dem_analysis_module = DEMAnalysisStage

        if self.settings["settings"].Has("working_directory"):
            working_dir = self.settings["settings"]["working_directory"].GetString()

            class DEMAnalysisStageWithWorkingDir(DEMAnalysisStage):
                @classmethod
                def GetMainPath(self):
                    return working_dir

            dem_analysis_module = DEMAnalysisStageWithWorkingDir

        return dem_analysis_module(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        super(DEMWrapper,self).SolveSolutionStep()

        for mp_name in self.settings["settings"]["move_mesh_model_part"].GetStringArray():
            DEMApplication.MoveMeshUtility().MoveDemMesh(self.model[mp_name].Nodes,True)
