import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem


class ControlModuleFemDemUtility(object):
    def __init__(self, Model, spheres_model_part):

        self.components_utility_list = []

        self.top_fem_model_part = Model["Structure.top_fem_surface"]
        self.top_dem_model_part = spheres_model_part.GetSubModelPart("bot_dem_surface")
        top_settings = KratosMultiphysics.Parameters( """
        {
            "variable_name": "DISPLACEMENT",
            "reaction_variable_name": "REACTION",
            "dem_force_variable_name": "TOTAL_FORCES",
            "target_stress_variable_name": "TARGET_STRESS",
            "reaction_stress_variable_name": "REACTION_STRESS",
            "imposed_direction" : 2,
            "target_stress_table_id" : 1,
            "initial_velocity" : -0.15,
            "limit_velocity" : -0.6,
            "velocity_factor" : 0.5,
            "compression_length" : 0.003,
            "young_modulus" : 7.0e9,
            "start_time" : 0.0,
            "face_area": 1.0
        }  """ )

        self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.top_fem_model_part,self.top_dem_model_part,top_settings))

        self.bot_fem_model_part = Model["Structure.bot_fem_surface"]
        self.bot_dem_model_part = spheres_model_part.GetSubModelPart("top_dem_surface")
        bot_settings = KratosMultiphysics.Parameters( """
        {
            "variable_name": "DISPLACEMENT",
            "reaction_variable_name": "REACTION",
            "dem_force_variable_name": "TOTAL_FORCES",
            "target_stress_variable_name": "TARGET_STRESS",
            "reaction_stress_variable_name": "REACTION_STRESS",
            "imposed_direction" : 2,
            "target_stress_table_id" : 2,
            "initial_velocity" : 0.15,
            "limit_velocity" : 0.6,
            "velocity_factor" : 0.5,
            "compression_length" : 0.003,
            "young_modulus" : 7.0e9,
            "start_time" : 0.0,
            "face_area": 1.0
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