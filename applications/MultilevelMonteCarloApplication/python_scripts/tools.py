from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics


class ParametersWrapper(object):
    """
    Class for handling the project parameters with different solver settings
    input:  project_parameters: Kratos parameters class before refinement
    """
    def __init__(self,project_parameters):
        self.project_parameters = project_parameters

    """
    function returning the model part name
    input:  self: an instance of the class
    output: model_part_name: main model part name
    """
    def GetModelPartName(self):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["model_part_name"].GetString()
        else:
            return project_parameters["solver_settings"]["model_part_name"].GetString()

    """
    function returning the domain size
    input:  self: an instance of the class
    output: domain_size: domain size
    """
    def GetDomainSize(self):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["domain_size"].GetInt()
        else:
            return project_parameters["solver_settings"]["domain_size"].GetInt()
    """
    function returning project parameters with modified input type of the model import settings
    input:  self:   an instance of the class
            string: string to be set as input type
    output: project_parameters: modified project parameters
    """
    def SetModelImportSettingsInputType(self,string):
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters
        else:
            project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters