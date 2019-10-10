from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

class ParametersWrapper(object):
    """
    input:  project_parameters: Kratos parameters class before refinement
    """
    def __init__(self,project_parameters):
        self.project_parameters = project_parameters

    def GetModelPartName(self):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["model_part_name"].GetString()
        else:
            return project_parameters["solver_settings"]["model_part_name"].GetString()

    def GetDomainSize(self):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["domain_size"].GetInt()
        else:
            return project_parameters["solver_settings"]["domain_size"].GetInt()

    def SetModelImportSettingsInputType(self,string):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters
        else:
            project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters