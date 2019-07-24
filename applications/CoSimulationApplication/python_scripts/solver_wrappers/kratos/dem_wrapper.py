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

        CurrentDEMAnalysisStage = DEMAnalysisStage(self.model, self.project_parameters)

        if "path_to_project" in self.settings["settings"].keys():
            def NewGetMainPath(self):
                return self.settings["settings"]["path_to_project"].GetString()
            import types
            CurrentDEMAnalysisStage.GetMainPath = types.MethodType(NewGetMainPath, CurrentDEMAnalysisStage)


        return CurrentDEMAnalysisStage

    def SolveSolutionStep(self):
        super(DEMWrapper,self).SolveSolutionStep()

        for mp_name in self.settings["settings"]["move_mesh_model_part"].GetStringArray():
            DEMApplication.MoveMeshUtility().MoveDemMesh(self.model[mp_name].Nodes,True)
