import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem


class MultiaxialControlModuleFEMDEMGeneralized2DUtility(object):
    def __init__(self, dem_model_part, Model):

        self.dem_model_part = dem_model_part
        self.fem_model_part = Model["Structure"]

        project_parameters_file_name = "SPParameters.json"

        with open(project_parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.parameters = project_parameters["multiaxial_control_module_fem_dem_generalized_2d_utility"]

        # Negative target_stress means compression.

        self.cm_utility = DemFem.MultiaxialControlModuleFEMDEMGeneralized2DUtilities(self.dem_model_part, 
                                                                                self.fem_model_part, 
                                                                                self.parameters)

    def ExecuteInitialize(self):
        self.cm_utility.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.cm_utility.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cm_utility.ExecuteFinalizeSolutionStep()