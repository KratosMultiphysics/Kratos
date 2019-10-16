import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem


class ControlModuleFemDemUtility(object):
    def __init__(self, Model, spheres_model_part, test_number):

        if not test_number:
            return

        self.components_utility_list = []

        compression_length = 0.00381
        if test_number == 1: # CTW16
            face_area = 0.008062
        elif test_number == 2: # CTW10
            face_area = 0.007601
        else: # Blind test
            compression_length = 0.009144
            face_area = 0.088343

        self.top_fem_model_part = Model["Structure.SurfacePressure3D_top_pressure"]
        #self.top_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZpos"]
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
            "target_stress_table_id" : 1,
            "initial_velocity" : 0.0,
            "limit_velocity" : -0.1,
            "velocity_factor" : 0.5,
            "young_modulus" : 7.0e9,
            "stress_increment_tolerance": 100.0,
            "update_stiffness": true,
            "start_time" : 0.0
        }  """ )

        top_settings.AddEmptyValue("compression_length")
        top_settings["compression_length"].SetDouble(compression_length)
        top_settings.AddEmptyValue("face_area")
        top_settings["face_area"].SetDouble(face_area)
        self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.top_fem_model_part, self.top_dem_model_part, top_settings))

        self.bot_fem_model_part = Model["Structure.SurfacePressure3D_bottom_pressure"]
        #self.bot_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZneg"]
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
            "target_stress_table_id" : 2,
            "initial_velocity" : 0.0,
            "limit_velocity" : 0.1,
            "velocity_factor" : 0.5,
            "young_modulus" : 7.0e9,
            "stress_increment_tolerance": 100.0,
            "update_stiffness": true,
            "start_time" : 0.0
        }  """ )

        bot_settings.AddEmptyValue("compression_length")
        bot_settings["compression_length"].SetDouble(compression_length)
        bot_settings.AddEmptyValue("face_area")
        bot_settings["face_area"].SetDouble(face_area)
        self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.bot_fem_model_part, self.bot_dem_model_part, bot_settings))

    def ExecuteInitialize(self):

        for component in self.components_utility_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_utility_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_utility_list:
            component.ExecuteFinalizeSolutionStep()