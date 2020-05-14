import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem


class MultiaxialControlModuleGeneralized2DUtility(object):
    def __init__(self, dem_model_part, dem_fem_boundary_model_part):

        self.dem_model_part = dem_model_part
        self.dem_fem_boundary_model_part = dem_fem_boundary_model_part

        project_parameters_file_name = "sp_2d_rigid_fem_parameters.json"

        with open(project_parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.parameters = project_parameters["multiaxial_control_module_generalized_2d_utility"]

        # Negative target_stress means compression.

        self.cm_utility = Dem.MultiaxialControlModuleGeneralized2DUtilities(self.dem_model_part, 
                                                                            self.dem_fem_boundary_model_part, 
                                                                            self.parameters)

    def ExecuteInitialize(self):
        self.cm_utility.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.cm_utility.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cm_utility.ExecuteFinalizeSolutionStep()