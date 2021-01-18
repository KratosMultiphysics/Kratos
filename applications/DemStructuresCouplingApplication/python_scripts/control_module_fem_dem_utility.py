import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem


class ControlModuleFemDemUtility():
    def __init__(self, Model, spheres_model_part, test_number):

        fem_main_model_part = Model["Structure"]
        self.dem_main_model_part = spheres_model_part

        self.components_utility_list = []

        if not test_number:
            return

        if test_number == 1: # CTW16
            compression_length = 0.00381
            # face_area = 0.008062
            alternate_axis_loading = False
            limit_velocity = -50.0
        elif test_number == 2: # CTW10
            compression_length = 0.00381
            # face_area = 0.007601
            alternate_axis_loading = False
            limit_velocity = -50.0
        else: # Blind test
            compression_length = 0.009144
            # face_area = 0.088343
            alternate_axis_loading = True
            limit_velocity = -15.0

        if fem_main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            self.fem_submodel_part = Model["Structure.Parts_Solid_part"]
            settings = KratosMultiphysics.Parameters( """
            {
                "target_stress_table_id" : 1,
                "initial_velocity" : 0.0,
                "velocity_factor" : 1.0,
                "compression_length" : 1.0,
                "young_modulus" : 7.0e9,
                "stress_increment_tolerance": 1.0e-3,
                "update_stiffness": true,
                "start_time" : 0.0,
                "stress_averaging_time": 1.0e-5
            }  """ )
            settings.AddEmptyValue("alternate_axis_loading")
            settings["alternate_axis_loading"].SetBool(alternate_axis_loading)
            settings.AddEmptyValue("limit_velocity")
            settings["limit_velocity"].SetDouble(limit_velocity)
            self.components_utility_list.append(DemFem.ControlModuleFemDem2DUtilities(self.fem_submodel_part, self.dem_main_model_part, settings))
        else:
            self.top_fem_model_part = Model["Structure.SurfacePressure3D_top_pressure"]
            #self.top_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZpos"]
            self.top_dem_model_part = self.dem_main_model_part.GetSubModelPart("topdem")
            top_settings = KratosMultiphysics.Parameters( """
            {
                "imposed_direction" : 2,
                "alternate_axis_loading": false,
                "target_stress_table_id" : 1,
                "initial_velocity" : 0.0,
                "limit_velocity" : -0.1,
                "velocity_factor" : 0.5,
                "young_modulus" : 7.0e9,
                "stress_increment_tolerance": 100.0,
                "update_stiffness": true,
                "start_time" : 0.0,
                "stress_averaging_time": 1.0e-5
            }  """ )

            top_settings.AddEmptyValue("compression_length")
            top_settings["compression_length"].SetDouble(compression_length)
            self.components_utility_list.append(DemFem.ControlModuleFemDemUtilities(self.top_fem_model_part, self.top_dem_model_part, top_settings))

            self.bot_fem_model_part = Model["Structure.SurfacePressure3D_bottom_pressure"]
            #self.bot_fem_model_part = Model["Structure.SurfacePressure3D_sigmaZneg"]
            self.bot_dem_model_part = self.dem_main_model_part.GetSubModelPart("botdem")
            bot_settings = KratosMultiphysics.Parameters( """
            {
                "imposed_direction" : 2,
                "alternate_axis_loading": false,
                "target_stress_table_id" : 2,
                "initial_velocity" : 0.0,
                "limit_velocity" : 0.1,
                "velocity_factor" : 0.5,
                "young_modulus" : 7.0e9,
                "stress_increment_tolerance": 100.0,
                "update_stiffness": true,
                "start_time" : 0.0,
                "stress_averaging_time": 1.0e-5
            }  """ )

            bot_settings.AddEmptyValue("compression_length")
            bot_settings["compression_length"].SetDouble(compression_length)
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