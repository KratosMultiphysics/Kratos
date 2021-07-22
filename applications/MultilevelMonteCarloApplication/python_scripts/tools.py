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
        """
        if (self.project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid"):
            self.project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]["input_type"].SetString(string)
        else:
            self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString(string)

    def GetMaterialsFilename(self):
        """
        Method returning materials filename.

        Input:
        - self: an instance of the class

        Output:
        - materials_filename: the materials filename
        """
        if self.project_parameters["solver_settings"].Has("material_import_settings"):
            return self.project_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
        else:
            return None


    def SetMaterialsFilename(self,string):
        """
        Method returning project parameters with modified materials filename.

        Input:
        - self: an instance of the class
        - string: string to be set as input type
        """
        if self.project_parameters["solver_settings"].Has("material_import_settings"):
            self.project_parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString(string)