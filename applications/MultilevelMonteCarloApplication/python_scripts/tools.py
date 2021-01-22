# Import Kratos
import KratosMultiphysics


class ParametersWrapper(object):
    """
    Class for handling the project parameters with different solver settings.
    This class is used by Monte Carlo and Multilevel Monte Carlo algorithms.

    Input:
    - project_parameters: Kratos parameters
    """
    def __init__(self,project_parameters):
        self.project_parameters = project_parameters

    def GetModelPartName(self):
        """
        Method returning the model part name.

        Input:
        - self: an instance of the class

        Output:
        - model_part_name: string containing the main model part name
        """
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["model_part_name"].GetString()
        else:
            return project_parameters["solver_settings"]["model_part_name"].GetString()

    def GetDomainSize(self):
        """
        Method returning the domain size of the problem

        Input:
        - self: an instance of the class

        Output:
        - domain_size: domain size of the problem
        """
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            return project_parameters["solver_settings"]["fluid_solver_settings"]["domain_size"].GetInt()
        else:
            return project_parameters["solver_settings"]["domain_size"].GetInt()

    def SetModelImportSettingsInputType(self,string):
        """
        Method returning project parameters with modified input type of the model_import_settings.

        Input:
        - self: an instance of the class
        - string: string to be set as input type

        Output:
        - project_parameters: Kratos parameters with modified model_import_settings
        """
        project_parameters = self.project_parameters
        if (project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters
        else:
            project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString(string)
            return project_parameters