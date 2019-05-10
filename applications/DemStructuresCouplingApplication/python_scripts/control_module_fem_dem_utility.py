import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem


class ControlModuleFemDemUtility(object):
    def __init__(self, Model, spheres_model_part):

        self.components_utility_list = []

        self.top_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZpos"]
        self.top_dem_model_part = spheres_model_part.GetSubModelPart("topdem")
        top_settings = KratosMultiphysics.Parameters( """
        {
            "variable_name": "DISPLACEMENT",
            "reaction_variable_name": "REACTION",
            "dem_force_variable_name": "TOTAL_FORCES",
            "target_stress_variable_name": "TARGET_STRESS",
            "reaction_stress_variable_name": "REACTION_STRESS",
            "loading_velocity_variable_name": "LOADING_VELOCITY",
            "imposed_direction" : 2,
            "target_stress_table_id" : 2,
            "initial_velocity" : 0.0,
            "limit_velocity" : -0.1,
            "velocity_factor" : 0.5,
            "compression_length" : 0.009144,
            "young_modulus" : 7.0e9,
            "start_time" : 0.0,
            "face_area": 0.08834
        }  """ )

        self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.top_fem_model_part,self.top_dem_model_part,top_settings))

        self.bot_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZneg"]
        self.bot_dem_model_part = spheres_model_part.GetSubModelPart("botdem")
        bot_settings = KratosMultiphysics.Parameters( """
        {
            "variable_name": "DISPLACEMENT",
            "reaction_variable_name": "REACTION",
            "dem_force_variable_name": "TOTAL_FORCES",
            "target_stress_variable_name": "TARGET_STRESS",
            "reaction_stress_variable_name": "REACTION_STRESS",
            "loading_velocity_variable_name": "LOADING_VELOCITY",
            "imposed_direction" : 2,
            "target_stress_table_id" : 1,
            "initial_velocity" : 0.0,
            "limit_velocity" : 0.1,
            "velocity_factor" : 0.5,
            "compression_length" : 0.009144,
            "young_modulus" : 7.0e9,
            "start_time" : 0.0,
            "face_area": 0.08834
        }  """ )

        self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.bot_fem_model_part,self.bot_dem_model_part,bot_settings))

    def ExecuteInitialize(self):

        for component in self.components_utility_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_utility_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_utility_list:
            component.ExecuteFinalizeSolutionStep()