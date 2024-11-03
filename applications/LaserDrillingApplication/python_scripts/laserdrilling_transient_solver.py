import os
import numpy as np
from scipy.interpolate import interp1d
import h5py
import time as timer
from sympy import *
from KratosMultiphysics.LaserDrillingApplication.fem_tools import SurfaceFEMProjector
from os import environ
environ['OMP_NUM_THREADS'] = '4'

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver

def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolver(model, custom_settings)

class LaserDrillingTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the settings in base class
        super().__init__(model, custom_settings)   

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        self.jump_between_pulses_counter += self.delta_time
        error_in_delta_time = abs(self.jump_between_pulses_counter - self.time_jump_between_pulses)
        numerical_error = 1e-16
        if self.jump_between_pulses_counter >= self.time_jump_between_pulses or error_in_delta_time < numerical_error:
            self.jump_between_pulses_counter = 0
            self.pulse_number += 1
            self.ResetTemperatureField()

            if self.print_debug_info:
                initial_system_energy = self.MonitorEnergy()
                print("\n\nEnergy before laser:", initial_system_energy, "J")

            self.ImposeLaserDeltaTemperature()

            if self.print_debug_info:
                print("self.Q:", self.Q, "J")
                expected_energy_after_laser = initial_system_energy + self.Q
                system_energy = self.MonitorEnergy()
                print("Expected energy after laser:", expected_energy_after_laser, "J")
                print("Actual energy after laser:", system_energy, "J")
                relative_error = 100.0 * (system_energy - expected_energy_after_laser) / expected_energy_after_laser
                print("Relative error in energy (%):", relative_error, "\n\n")

            self.RemoveElementsByAblation()
            self.AdjustTemperatureFieldAfterAblation()
            self.ResidualHeatStage()

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "time_integration_method" : "implicit",
            "transient_parameters" : {
                "dynamic_tau": 1.0,
                "theta"    : 0.5
            },
            "ambient_temperature" : 0.0
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AllocateKratosMemory(self):
        # Set element counter variable to zero
        for elem in self.main_model_part.Elements:
            elem.SetValue(LaserDrillingApplication.THERMAL_COUNTER, 0)
            elem.SetValue(KratosMultiphysics.TEMPERATURE, 0.0)
            elem.SetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE, 0.0)
            elem.SetValue(LaserDrillingApplication.THERMAL_DECOMPOSITION, 0.01)
            elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, 0.0)
            elem.SetValue(LaserDrillingApplication.ELEMENTAL_VOLUME, 0.0)
            elem.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, 0.0)
            elem.SetValue(LaserDrillingApplication.ENERGY_PER_VOLUME, 0.0)
            elem.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, 0.0)
            elem.Set(KratosMultiphysics.ACTIVE, True)
            elem.SetValue(LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME, 0.0)

        for node in self.main_model_part.Nodes:
            node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 0.0)
            node.SetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE, 0.0)
            node.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, 0.0)
            node.SetValue(LaserDrillingApplication.ENERGY_PER_VOLUME, 0.0)
            node.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, 0.0)
            node.SetValue(LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME, 0.0)

    def ComputeSpotDiameter(self):
        spot_diameter = self.beam_waist_diameter * sqrt(1.0 + ((self.focus_z_offset + self.z_ast_max) / self.rayleigh_length)**2)
        return spot_diameter

    def ComputePeakFluence(self):
        return 2.0 * self.Q / (np.pi * self.omega_0**2) # J/mm2

    def ComputeMaximumAblationRadius(self):
        import math
        return self.omega_0 * math.sqrt(0.5 * math.log(self.F_p / (self.delta_pen * self.q_ast)))

    def SetParameters(self):

        self.some_elements_are_above_the_evap_temp = False
        self.jump_between_pulses_counter = 0
        self.pulse_number = 0
        self.print_hdf5_and_gnuplot_files = False

        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        self.material_settings = materials["properties"][0]["Material"]

        with open("ProjectParameters.json", 'r') as project_parameters_file:
            self.project_parameters = KratosMultiphysics.Parameters(project_parameters_file.read())

        if not self.project_parameters["problem_data"].Has("average_laser_power"):
            self.average_laser_power = 18
        else:
            self.average_laser_power = self.project_parameters["problem_data"]["average_laser_power"].GetDouble()

        if not self.project_parameters["problem_data"].Has("pulse_frequency"):
            self.pulse_frequency = 2e5
        else:
            self.pulse_frequency = self.project_parameters["problem_data"]["pulse_frequency"].GetDouble()

        if not self.project_parameters["problem_data"].Has("beam_waist_diameter"):
            self.beam_waist_diameter = 0.0179
        else:
            self.beam_waist_diameter = self.project_parameters["problem_data"]["beam_waist_diameter"].GetDouble()

        if not self.project_parameters["problem_data"].Has("Rayleigh_length"):
            self.rayleigh_length = 0.409
        else:
            self.rayleigh_length = self.project_parameters["problem_data"]["Rayleigh_length"].GetDouble()

        if not self.project_parameters["problem_data"].Has("focus_Z_offset"):
            self.focus_z_offset = 0.4
        else:
            self.focus_z_offset = self.project_parameters["problem_data"]["focus_Z_offset"].GetDouble()

        self.z_ast_max = 0.0

        #self.spot_diameter = self.ComputeSpotDiameter()

        if not self.project_parameters["problem_data"].Has("vaporisation_temperature"):
            self.T_e = 1000.0
        else:
            self.T_e = self.project_parameters["problem_data"]["vaporisation_temperature"].GetDouble()

        if not self.project_parameters["problem_data"].Has("compute_vaporisation"):
            self.compute_vaporisation = False
        else:
            self.compute_vaporisation = self.project_parameters["problem_data"]["compute_vaporisation"].GetBool()

        if not self.material_settings["Variables"].Has("IONIZATION_ALPHA"):
            self.ionization_alpha = 1.0 #0.95
        else:
            self.ionization_alpha = 1.0 #self.material_settings['Variables']['IONIZATION_ALPHA'].GetDouble()

        if not self.material_settings["Variables"].Has("PENETRATION_DEPTH"):
            self.l_s = 0.002148 # mm.
        else:
            self.l_s = self.material_settings['Variables']['PENETRATION_DEPTH'].GetDouble()

        if not self.material_settings["Variables"].Has("ABLATION_THRESHOLD"):
            self.F_th = 0.010667 # J/mm2
        else:
            self.F_th = self.material_settings['Variables']['ABLATION_THRESHOLD'].GetDouble()

        if not self.material_settings["Variables"].Has("THERMAL_DEPTH"):
            self.l_th = 0.0007 # mm.
        else:
            self.l_th = self.material_settings['Variables']['THERMAL_DEPTH'].GetDouble()

        if not self.material_settings["Variables"].Has("ENTHALPY"):
            self.H_ev = 4e5 # J/Kg. Value found on the internet for a given epoxy resin.
        else:
            self.H_ev = self.material_settings['Variables']['ENTHALPY'].GetDouble()

        if not self.project_parameters["problem_data"].Has("mesh_size"):
            self.mesh_size = "coarse"
        else:
            self.mesh_size = self.project_parameters["problem_data"]["mesh_size"].GetString()

        if not self.project_parameters["problem_data"].Has("mesh_type"):
            self.mesh_type = "unstructured"
        else:
            self.mesh_type = self.project_parameters["problem_data"]["mesh_type"].GetString()

        if not self.project_parameters["problem_data"].Has("print_hole_geometry_files"):
            self.print_hole_geometry_files = False
        else:
            self.print_hole_geometry_files = self.project_parameters["problem_data"]["print_hole_geometry_files"].GetBool()

        self.Q = self.average_laser_power / self.pulse_frequency
        self.time_jump_between_pulses = 1.0 / self.pulse_frequency
        self.cp = self.material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = self.material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = self.material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        self.kappa = self.conductivity / (self.rho * self.cp)
        self.ablation_energy_fraction = self.ionization_alpha
        self.evaporation_energy_fraction = 1.0 - self.ionization_alpha

        ## 2024 Woodfield - Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue
        ## Laser data
        self.omega_0 = 0.5 * self.ComputeSpotDiameter() #self.R_far # mm
        import numpy as np
        y_limit = 2.0 * self.omega_0
        self.hole_theoretical_Y_coords = np.linspace(0.0, float(y_limit), 101)
        self.hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_Y_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)
        self.F_p = self.ComputePeakFluence()

        ## Material calibration using experiments
        if not self.material_settings["Variables"].Has("OPTICAL_PENETRATION_DEPTH"):
            self.delta_pen = 5e-4 # mm
        else:
            self.delta_pen = self.material_settings['Variables']['OPTICAL_PENETRATION_DEPTH'].GetDouble()

        if not self.material_settings["Variables"].Has("ENERGY_PER_VOLUME_THRESHOLD"):
            self.q_ast = 10.0 # J/mm3
        else:
            self.q_ast = self.material_settings['Variables']['ENERGY_PER_VOLUME_THRESHOLD'].GetDouble()

        if not self.material_settings["Variables"].Has("REFRACTIVE_INDEX"):
            self.refractive_index_n = 1.5
        else:
            self.refractive_index_n = self.material_settings['Variables']['REFRACTIVE_INDEX'].GetDouble()

        if not self.project_parameters["problem_data"].Has("consider_material_refraction"):
            self.consider_material_refraction = False
        else:
            self.consider_material_refraction = self.project_parameters["problem_data"]["consider_material_refraction"].GetBool()

        if self.material_settings['compute_optical_penetration_depth_using_refractive_index'].GetBool():
            self.ComputeOpticalPenetrationDepth()
            self.delta_pen = self.l_s

        if self.material_settings['compute_energy_per_unit_volume_threshold_using_enthalpy_and_ionization'].GetBool():
            self.ionizarion_energy_per_volume_threshold = self.ComputeIonizationEnergyPerUnitVolumeThreshold()
            self.use_enthalpy_and_ionization = True
        else:
            self.use_enthalpy_and_ionization = False

        self.decomposed_nodes_coords_filename = "hole_coords_q_ast=" + str(self.q_ast) + "_delta_pen=" + str(self.delta_pen) + "_" + self.mesh_type + "_" + self.mesh_size + ".txt"

        self.r_ast_max = self.ComputeMaximumAblationRadius()

        if not self.project_parameters["problem_data"].Has("adjust_T_field_after_ablation"):
            self.adjust_T_field_after_ablation = False
        else:
            self.adjust_T_field_after_ablation = self.project_parameters["problem_data"]["adjust_T_field_after_ablation"].GetBool()

        if not self.project_parameters["problem_data"].Has("reference_T_after_laser"):
            self.reference_T_after_laser = 298.15
        else:
            self.reference_T_after_laser = self.project_parameters["problem_data"]["reference_T_after_laser"].GetDouble()

        if not self.project_parameters["problem_data"].Has("print_debug_info"):
            self.print_debug_info = False
        else:
            self.print_debug_info = self.project_parameters["problem_data"]["print_debug_info"].GetBool()

        self.analytical_ablated_volume_in_n_pulses = 0.0

        for properties, i in enumerate(materials["properties"]):
            full_material_part_name_i = i["model_part_name"].GetString()
            prefix = 'ThermalModelPart.'
            material_part_name_i = full_material_part_name_i.removeprefix(prefix)

            material_part_i = self.main_model_part.GetSubModelPart(material_part_name_i)
            material_settings_i = i["Material"]
            thermal_energy_per_volume_i = material_settings_i['Variables']['ENERGY_PER_VOLUME_THRESHOLD'].GetDouble()
            for elem in material_part_i.Elements:
                elem.SetValue(LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME, thermal_energy_per_volume_i)

        #self.sigma = 0.5 * self.R_far
        #self.K = 1 / (2 * self.sigma**2)
        #import math
        #self.C = self.ablation_energy_fraction * self.Q * self.K / (np.pi * (1 - math.exp(-self.K * self.R_far**2)))
        #self.irradiated_surface_area = np.pi * self.R_far**2

        #if self.ablation_energy_fraction:
            # Find F_th_fraction multiplying F_th so Radius_th = R_far
            # This gives a F_th_fraction of 0.313, a little too flat (maximum not captured)
            # F_th_fraction = self.C * np.exp(-self.K * self.R_far**2) / (self.ablation_energy_fraction * self.Q / A)

            ####################################
            # Calibrated for:                  #
            # residual heat fraction   = 0.05  #
            # Vaporisation temperature = 1000K #
            # Power                    = 3W    #
            ####################################

            #self.F_th = 0.009667 # J/mm2 #F_th_fraction * self.ablation_energy_fraction * self.Q / self.irradiated_surface_area
            #self.radius_th = math.sqrt(math.log(self.C / self.F_th) / self.K)

            #self.l_s = self.PenetrationDepthEstimation() #0.001 # 5 pulses, with evaporation

        #self.l_s = 0.002048
        #self.l_th = 0.33 # micrometers # Thermal depth, assumed

        l_th_in_meters = self.l_th * 1e-3
        kappa_in_square_meters = self.kappa * 1e-6
        self.thermal_penetration_time = l_th_in_meters**2 / kappa_in_square_meters

        # Finite Elements
        #self.n_surface_elements = 10 # number of elements
        #self.sparse_option = True
        #self.surface_nodes_Y_values = np.linspace(0.0, self.R_far, self.n_surface_elements + 1)
        #self.surface_element_size = self.R_far / self.n_surface_elements

        # Debug
        self.plot_progressive_hole_figures = False

    def ComputeMaximumDepth(self):
        maximum_depth = self.list_of_ablated_nodes_coords_X[0]
        return maximum_depth

    def ComputeOpticalPenetrationDepth(self):
        light_lambda = 550e-6 # mm, light wavelength
        epoxy_n = self.refractive_index_n
        n = epoxy_n
        A = 4.0 * n / ((n + 1)**2 + n**2)
        self.l_s = 0.25 * light_lambda * A / np.pi
        return self.l_s

    def UpdateLaserRelatedParameters(self):
        self.z_ast_max = self.ComputeMaximumDepth()
        self.omega_0 = 0.5 * self.ComputeSpotDiameter()
        self.F_p = self.ComputePeakFluence()
        self.r_ast_max = self.ComputeMaximumAblationRadius()

    def SetUpResultsFiles(self):
        self.SetUpGNUPlotFiles()
        self.SetUpHDF5Files()

    def SetUpGNUPlotFiles(self):
        if os.path.exists("temperature_alpha.txt"):
            os.remove("temperature_alpha.txt")
        self.temperature_alpha_file = open("temperature_alpha.txt", "a")
        if os.path.exists("time_alpha.txt"):
            os.remove("time_alpha.txt")
        self.time_alpha_file = open("time_alpha.txt", "a")
        if os.path.exists("decomposed_volume_evolution.txt"):
            os.remove("decomposed_volume_evolution.txt")
        self.decomposed_volume_file = open("decomposed_volume_evolution.txt", "a")

    def SetUpHDF5Files(self):
        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2
        self.near_field_nodes = [node for node in self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])        

    def AdjustTemperatureFieldAfterAblation(self):
        if not self.adjust_T_field_after_ablation:
            return
        for node in self.main_model_part.Nodes:
            old_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            new_temperature = old_temperature + self.T0 - self.reference_T_after_laser
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temperature)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temperature)
        for elem in self.main_model_part.Elements:
            elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)

    def ImposeLaserDeltaTemperature(self):
        if not self.consider_material_refraction:
            self.ImposeTemperatureIncreaseDueToLaser()
        else:
            self.ImposeTemperatureIncreaseDueToLaserWithRefraction()

    def Initialize(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).Initialize()
        self.starting_time = timer.time()

        self.AllocateKratosMemory()

        self.SetParameters()

        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []
        self.list_of_lists_of_decomposed_nodes_X = []
        self.list_of_lists_of_decomposed_nodes_Y = []
        self.full_list_of_ablated_nodes_coords_X = []
        self.full_list_of_ablated_nodes_coords_Y = []

        if self.print_hdf5_and_gnuplot_files:
            self.SetUpResultsFiles()

        # TODO: Initial condition, ambient temperature. Do this using GUI!
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
            #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)

        self.IdentifyInitialSurfaceNodes()

        self.ResidualHeatStage()

        if self.print_hdf5_and_gnuplot_files:
            self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def CountActiveElements(self):
        count = 0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                count += 1
        return count

    def IdentifyInitialSurfaceNodes(self):
        list_of_decomposed_nodes_ids = []
        list_of_decomposed_elements_ids = []            
        list_of_decomposed_nodes_coords_X = []
        list_of_decomposed_nodes_coords_Y = []

        for elem in self.main_model_part.Elements:
            for node in elem.GetNodes():
                # TODO: extremely ad-hoc!
                if node.X < 0.0000001 and node.Y <= self.r_ast_max:
                    list_of_decomposed_nodes_ids.append(node.Id)
                    list_of_decomposed_elements_ids.append(elem.Id)    
                    list_of_decomposed_nodes_coords_X.append(node.X)
                    list_of_decomposed_nodes_coords_Y.append(node.Y)                

        self.list_of_decomposed_nodes_ids = np.array(list(set(list_of_decomposed_nodes_ids)))
        self.list_of_decomposed_elements_ids = np.array(list(set(list_of_decomposed_elements_ids)))
        self.list_of_decomposed_nodes_coords_X = np.array(list_of_decomposed_nodes_coords_X)
        self.list_of_decomposed_nodes_coords_Y = np.array(list_of_decomposed_nodes_coords_Y)

        if not self.main_model_part.HasSubModelPart("BoundaryPart"):
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')

        self.boundary_part.AddElements(self.list_of_decomposed_elements_ids)
        self.boundary_part.AddNodes(self.list_of_decomposed_nodes_ids)

    def StorePreEvaporationTemperature(self):
        for node in self.main_model_part.Nodes:
            pre_evaporation_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            node.SetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE, pre_evaporation_temperature)
            #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)
        for elem in self.main_model_part.Elements:
            elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)
            pre_evaporation_temperature = elem.GetValue(KratosMultiphysics.TEMPERATURE)
            elem.SetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE, pre_evaporation_temperature)

    def ResetTemperatureField(self):
        if self.adjust_T_field_after_ablation:
            reference_temp = self.reference_T_after_laser
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, reference_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, reference_temp)

    def ResidualHeatStage(self):
        if self.evaporation_energy_fraction:
            self.projector = SurfaceFEMProjector(self.n_surface_elements, self.R_far, self.sparse_option) #, delta_coefficients)
            self.q_interp = self.projector.InterpolateFunctionAndNormalize(self.EnergyPerUnitArea1D) #, 1.0)
            if self.compute_vaporisation:
                self.first_evaporation_stage_done = False
                self.max_vaporisation_layers = 50
                # if self.pulse_number == 1:
                #     self.max_vaporisation_layers = 1
                self.vaporisation_layer_number = 1
                self.some_elements_are_above_the_evap_temp = True
                print("Removing elements by evaporation...")
                print("Pulse number:", self.pulse_number)
                self.last_evaporation_layer_applied = False
                self.StorePreEvaporationTemperature()
                self.ablation_only_legend_added = False

                while self.some_elements_are_above_the_evap_temp:
                    if self.vaporisation_layer_number > self.max_vaporisation_layers:
                        if not self.max_vaporisation_layers:
                            self.ImposeTemperatureIncreaseDueTo1DConduction()
                        print("******************************MAXIMUM ITERATIONS EXCEEDED!!!")
                        break
                    self.ImposeTemperatureIncreaseDueTo1DConduction()
                    self.RemoveElementsByEvaporation()

                print("Done!")
            else:
                self.ImposeTemperatureIncreaseDueTo1DConduction()

    def RemoveElementsByEvaporation(self):
        evap_elements_centers_Y = []
        evap_elements_volumes   = []
        delta_temp_elements     = []
        uncapped_delta_temp_elements     = []
        self.some_elements_are_above_the_evap_temp = False
        for elem in self.boundary_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                temp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)
                element_temperature = temp[0]
                if element_temperature > self.T_e:
                    pre_evap_temp = elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE)
                    delta_temp = self.T_e - pre_evap_temp
                    uncapped_delta_temp = element_temperature - pre_evap_temp
                    delta_temp_elements.append(delta_temp)
                    uncapped_delta_temp_elements.append(uncapped_delta_temp)
                    self.some_elements_are_above_the_evap_temp = True
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    vol = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
                    elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, vol[0])
                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
                    Y_centroid = elem.GetGeometry().Center().Y
                    evap_elements_centers_Y.append(Y_centroid)
                    evap_elements_volumes.append(element_volume)
                    for node in elem.GetNodes():
                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)

        number_of_problematic_elements = 0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                number_of_decomposed_nodes = 0
                temp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)
                element_temperature = temp[0]
                for node in elem.GetNodes():
                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                        number_of_decomposed_nodes +=1
                if number_of_decomposed_nodes == 3:
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
                    Y_centroid = elem.GetGeometry().Center().Y
                    pre_evap_temp = elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE)
                    delta_temp = self.T_e - pre_evap_temp
                    #delta_temp_elements.append(delta_temp)
                    uncapped_delta_temp = element_temperature - pre_evap_temp
                    #uncapped_delta_temp_elements.append(uncapped_delta_temp)
                    #evap_elements_centers_Y.append(Y_centroid)
                    #evap_elements_volumes.append(element_volume)                    
                    number_of_problematic_elements += 1
        print("Number_of_problematic_elements:", number_of_problematic_elements)
        evap_elements_centers_Y = np.array(evap_elements_centers_Y)
        evap_elements_volumes   = np.array(evap_elements_volumes)
        self.evap_elements_centers_Y = evap_elements_centers_Y[evap_elements_centers_Y.argsort()]
        self.evap_elements_volumes = evap_elements_volumes[evap_elements_centers_Y.argsort()]
        #print("self.evap_elements_centers_Y:", self.evap_elements_centers_Y)
        #print("self.evap_elements_volumes:", self.evap_elements_volumes)

        delta_temp_elements     = np.array(delta_temp_elements)
        self.delta_temp_elements = delta_temp_elements[evap_elements_centers_Y.argsort()]
        uncapped_delta_temp_elements = np.array(uncapped_delta_temp_elements)
        self.uncapped_delta_temp_elements = uncapped_delta_temp_elements[evap_elements_centers_Y.argsort()]

        # Total enthalpy: Energy consumed for the material to vaporize plus the energy for heating up the material to the vaporization temperature.
        # Equation (8) in Wang, 2019. 'Thermal effect of femtosecond laser polystyrene processing'
        self.evap_elements_enthalpies = self.evap_elements_volumes * self.rho * (self.H_ev + self.cp * self.delta_temp_elements)

        #print("evap_elements_volumes:", evap_elements_volumes)
        #print("delta_temp_elements:", delta_temp_elements)
        #print("self.evap_elements_enthalpies:", self.evap_elements_enthalpies)

        self.support_elements = [[] for i in range(self.n_surface_elements+1)]

        if not self.some_elements_are_above_the_evap_temp:
            return

        import matplotlib.pyplot as plt
        figure, axis = plt.subplots(2, 2, figsize=(15,8))
        label_size = 13
        numbers_size = 11
        axis[0][0].grid()
        axis[0][0].plot(self.evap_elements_centers_Y, self.evap_elements_volumes, color='red', marker='+')
        axis[0][0].set_xlabel('radius (mm)')
        axis[0][0].set_ylabel("Volumes (mm3)", fontsize = label_size)
        axis[1][0].grid()
        axis[1][0].plot(self.evap_elements_centers_Y, self.delta_temp_elements, color='blue', marker='o')
        #axis[1][0].set_ylim(bottom=0, top=None)
        axis[1][0].set_xlabel('radius (mm)')
        axis[1][0].set_ylabel("Delta temps (K)", fontsize = label_size)
        axis[0][1].grid()
        axis[0][1].plot(self.evap_elements_centers_Y, self.evap_elements_enthalpies, color='black', marker='x')
        axis[0][1].set_xlabel('radius (mm)')
        axis[0][1].set_ylabel("Enthalpies (J)", fontsize = label_size)
        axis[1][1].grid()
        axis[1][1].plot(self.evap_elements_centers_Y, self.uncapped_delta_temp_elements, color='green', marker='*')
        #axis[1][1].set_ylim(bottom=0, top=None)
        axis[1][1].set_xlabel('radius (mm)')
        axis[1][1].set_ylabel("Uncapped delta temps (K)", fontsize = label_size)
        #figure.show()

        print("\nEvaporating layer number", self.vaporisation_layer_number, "...")
        self.vaporisation_layer_number += 1

        if not self.sparse_option:
            self.projector.FillUpMassMatrix()
        else:
            self.projector.FillUpSparseMassMatrix()
        self.projector.AssignDeltasToTestFunctionSupports(self.evap_elements_centers_Y, self.support_elements)
        self.projector.FillUpDeltasRHS(self.evap_elements_centers_Y, self.support_elements, self.evap_elements_enthalpies)
        self.u = self.projector.Project()

        total_energy = self.projector.CalculateEnergyOfFEMFunction(self.u)

        print('Total energy expected =', sum(self.evap_elements_enthalpies))
        print('Total energy calculated =', total_energy)

        import matplotlib.pyplot as plt
        _, axis = plt.subplots(1, 2, figsize=(20,5))
        axis[0].grid()
        axis[0].plot(self.projector.X, self.q_interp, color='red', marker='+')
        axis[0].plot(self.projector.X, self.u, color='blue', marker='o')

        self.q_interp -= self.u

        # TODO: rethink this!
        for i, q in enumerate(self.q_interp):
            if q < 0.0:
                self.q_interp[i] = 0.0

        axis[0].tick_params(axis='both', which='major', labelsize = numbers_size)
        axis[0].tick_params(axis='both', which='minor', labelsize = numbers_size)
        plt.subplots_adjust(left=0.1, bottom=None, right=0.9, top=None, wspace=0.22, hspace=None)
        axis[0].plot(self.projector.X, self.q_interp, color='black', marker='*')
        axis[0].set_ylim(bottom=0.0, top=0.0036)
        axis[0].legend(["fluence (interpolated)", "fluence (lost)", "fluence (remaining)"], loc="upper right", fontsize = label_size)
        axis[0].set_xlabel('radius (mm)', fontsize = label_size)
        axis[0].set_ylabel("Energies (J/mm2)", fontsize = label_size) 

        self.AddDecomposedNodesToSurfaceList()
        self.first_evaporation_stage_done = True

        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y

        if not self.ablation_only_legend_added:
            if self.ablation_energy_fraction:
                self.full_list_of_ablated_nodes_coords_X.append(self.list_of_ablated_nodes_coords_X)
                self.full_list_of_ablated_nodes_coords_Y.append(self.list_of_ablated_nodes_coords_Y)
            else:
                self.full_list_of_ablated_nodes_coords_X.append(0.0*Y)
                self.full_list_of_ablated_nodes_coords_Y.append(Y)
            self.ablation_only_legend_added = True

        self.list_of_lists_of_decomposed_nodes_X.append(X)
        self.list_of_lists_of_decomposed_nodes_Y.append(Y)

        axis[1].grid()
        list_of_legends = []
        p = 1
        axis[1].tick_params(axis='both', which='major', labelsize = numbers_size)
        axis[1].tick_params(axis='both', which='minor', labelsize = numbers_size)
        for list_X, list_Y in zip(self.full_list_of_ablated_nodes_coords_X, self.full_list_of_ablated_nodes_coords_Y):
            axis[1].plot(list_Y, -list_X, color='black')
            axis[1].set_ylim(bottom=-0.0023, top=None)
            converted_num = "Ablation number #" + str(p)
            list_of_legends.append(converted_num)
            p += 1
        i = 1
        for list_X, list_Y in zip(self.list_of_lists_of_decomposed_nodes_X, self.list_of_lists_of_decomposed_nodes_Y):
            axis[1].plot(list_Y, -list_X)
            axis[1].set_ylim(bottom=-0.0023, top=0)
            converted_num = "Evap. layer #" + str(i)
            list_of_legends.append(converted_num)
            i += 1
        axis[1].legend(list_of_legends, loc="lower right", fontsize = 11)
        axis[1].set_xlabel('radius (mm)', fontsize = label_size)
        axis[1].set_ylabel("Ablation plus evaporation (mm)", fontsize = label_size) 

        if self.plot_progressive_hole_figures:
            plt.show()

        if self.vaporisation_layer_number <= self.max_vaporisation_layers:
            self.RetrievePreEvaporationTemperatureState()

    def RetrievePreEvaporationTemperatureState(self):
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, node.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, node.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE))
        for elem in self.main_model_part.Elements:
            elem.SetValue(KratosMultiphysics.TEMPERATURE, elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE))

    def ImposeTemperatureIncreaseDueTo1DConduction(self):
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        self.minimum_characteristic_Z =  1e6
        self.maximum_characteristic_Z = -1e6
        for node in self.main_model_part.Nodes:
            radius = node.Y
            """ if not self.ablation_energy_fraction:
                distance_to_surface = 0.0
            else: """
            F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
            distance_to_surface = F(radius)
            z = node.X - distance_to_surface
            if radius <= self.r_ast_max:
                delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
                old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temp = old_temp + delta_temp
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)
        #print("\nResidual heat fraction:", self.ionization_alpha)
        #print("\nMaximum characteristic depth:", self.maximum_characteristic_Z)
        problem_characteristic_time_minimum_depth = self.minimum_characteristic_Z**2 / self.kappa
        problem_characteristic_time_maximum_depth = self.maximum_characteristic_Z**2 / self.kappa
        minimum_time_step_for_minimum_depth = 0.1 * problem_characteristic_time_minimum_depth
        minimum_time_step_for_maximum_depth = 0.1 * problem_characteristic_time_maximum_depth
        #print("\nThermal problem characteristic time for maximum depth:", problem_characteristic_time_maximum_depth)
        #print("\nNecessary time step for maximum depth:", minimum_time_step_for_maximum_depth, '\n')

    def EnergyPerUnitArea1D(self, radius):
        C = (1 - self.ablation_energy_fraction) * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
        q =  C * np.exp(-self.K * radius**2)
        return q

    def InitialThermalConductionTime(self, radius):
        # This function returns the characteristic time required for the initial heat distribution from the surface to the interior
        # 4.0 (2**2) in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        # C = 4.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * (self.T_e - self.T0)**2)
        # t = C * self.EnergyPerUnitArea1D(radius)**2
        return self.thermal_penetration_time

    def TemperatureVariationInZDueToLaser1D(self, radius, z):
        #q = self.EnergyPerUnitArea1D(radius)
        q = self.projector.EvaluateFEMFunction(self.q_interp, radius)
        t_penetration = self.InitialThermalConductionTime(radius)
        # 2.0 in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        C = 2.0 * q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_penetration))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_penetration))
        characteristic_Z = np.sqrt(4.0 * self.kappa * t_penetration)
        if characteristic_Z <= self.minimum_characteristic_Z:
            self.minimum_characteristic_Z = characteristic_Z
        if characteristic_Z >= self.maximum_characteristic_Z:
            self.maximum_characteristic_Z = characteristic_Z
        return delta_temp

    def SolveSolutionStep(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).SolveSolutionStep()
        if self.print_hdf5_and_gnuplot_files:
            self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
            self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def ComputePulseHoleAndAddToTotalHole(self):
        for i, Y_coord in enumerate(self.hole_theoretical_Y_coords):
            self.hole_theoretical_X_coords[i] += self.EvaporationDepth(Y_coord)
        if self.pulse_number == 1:
            import copy
            self.one_pulse_hole_theoretical_X_coords = copy.deepcopy(self.hole_theoretical_X_coords)
            self.one_pulse_hole_theoretical_Y_coords = self.hole_theoretical_Y_coords

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.print_hdf5_and_gnuplot_files:
            decomp_vol = self.MonitorDecomposedVolume()
            decomp_vol *= 1e9 # To convert mm3 into um3
            current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            self.decomposed_volume_file.write(str(current_time) + " " + str(decomp_vol) + "\n")
            for elem in self.main_model_part.Elements:
                if elem.Id == self.element_id_to_study:
                    temperature = elem.GetValue(KratosMultiphysics.TEMPERATURE)
                    thermal_decomposition = elem.GetValue(LaserDrillingApplication.THERMAL_DECOMPOSITION)
                    self.temperature_alpha_file.write(str(temperature) + " " + str(thermal_decomposition) + "\n")
                    self.time_alpha_file.write(str(current_time) + " " + str(thermal_decomposition) + "\n")
                    break
        self.UpdateLaserRelatedParameters()
        self.ComputePulseHoleAndAddToTotalHole()

    def Finalize(self):
        super().Finalize()
        elapsed_time = timer.time() - self.starting_time
        print("\nElapsed_time:", elapsed_time, '\n')
        if self.print_hdf5_and_gnuplot_files:
            self.decomposed_volume_file.close()
            self.temperature_alpha_file.close()
            self.time_alpha_file.close()
            self.PrintDecomposedVolumeEvolution()

    def PrintDecomposedVolumeEvolution(self):
        import matplotlib.pyplot as plt
        file = open(os.path.expanduser("decomposed_volume_evolution.txt"))
        lines = file.readlines()
        x, y = [], []
        for line in lines:
            x.append(line.split()[0])
            y.append(line.split()[1])
        file.close()
        plt.plot(x,y)
        plt.show()

    def MonitorDecomposedVolume(self):
        decomposed_volume = 0.0
        for elem in self.main_model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo)
            decomposed_volume += out[0]
        return decomposed_volume

    def MonitorEnergy(self):
        energy = 0.0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
                out = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.THERMAL_ENERGY, self.main_model_part.ProcessInfo)
                energy += out[0]
        return energy

    def RemoveElementsByAblation(self):
        if self.ablation_energy_fraction:
            X = self.list_of_decomposed_nodes_coords_X
            Y = self.list_of_decomposed_nodes_coords_Y

            for elem in self.main_model_part.Elements:
                X_centroid = elem.GetGeometry().Center().X
                Y_centroid = elem.GetGeometry().Center().Y
                # if self.pulse_number == 1:
                #     X_interp = 0
                # else:
                F = interp1d(Y, X, bounds_error=False)
                X_interp = F(Y_centroid)
                DeltaX = X_centroid - X_interp
                d_ev = self.EvaporationDepth(Y_centroid)
                if DeltaX <= d_ev: # and Y_centroid <= self.radius_th:
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    for node in elem.GetNodes():
                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)
            for elem in self.main_model_part.Elements:
                if elem.Is(KratosMultiphysics.ACTIVE):
                    number_of_decomposed_nodes = 0
                    for node in elem.GetNodes():
                        if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                            number_of_decomposed_nodes += 1
                    if number_of_decomposed_nodes == 3:
                        elem.Set(KratosMultiphysics.ACTIVE, False)
            self.AddDecomposedNodesToSurfaceList()
            self.list_of_ablated_nodes_coords_X = self.list_of_decomposed_nodes_coords_X
            self.list_of_ablated_nodes_coords_Y = self.list_of_decomposed_nodes_coords_Y
            '''print('\nR_far:', self.R_far)
            print('\nRadius_th:', self.radius_th)
            print("\nDecomposed volume:", self.MonitorDecomposedVolume())'''

    def sortSecond(self, val):
        return val[1]

    def AddDecomposedNodesToSurfaceList(self):
        list_of_decomposed_nodes_ids = []
        list_of_decomposed_elements_ids = []
        number_of_boundary_elements = 0
        for elem in self.main_model_part.Elements:
            first_decomposed_node_found = False
            if elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                        list_of_decomposed_nodes_ids.append(node.Id)
                        if not first_decomposed_node_found:
                            number_of_boundary_elements += 1
                            list_of_decomposed_elements_ids.append(elem.Id)
                            first_decomposed_node_found = True
        self.list_of_decomposed_nodes_ids = np.array(list(set(list_of_decomposed_nodes_ids)))
        self.list_of_decomposed_elements_ids = np.array(list_of_decomposed_elements_ids)

        if not self.main_model_part.HasSubModelPart("BoundaryPart"):
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')
        else:
            self.main_model_part.RemoveSubModelPart(self.boundary_part)
            self.main_model_part.CreateSubModelPart('BoundaryPart')
            self.boundary_part = self.main_model_part.GetSubModelPart('BoundaryPart')
        self.boundary_part.AddElements(self.list_of_decomposed_elements_ids)
        self.boundary_part.AddNodes(self.list_of_decomposed_nodes_ids)

        list_of_decomposed_nodes_coords = []
        for node in self.main_model_part.Nodes:
            if node.Id in self.list_of_decomposed_nodes_ids:
                X = node.X
                Y = node.Y
                coords = [X, Y]
                list_of_decomposed_nodes_coords.append(coords)
        list_of_decomposed_nodes_coords.sort(key=self.sortSecond)
        self.list_of_decomposed_nodes_coords_X = np.array([coord[0] for coord in list_of_decomposed_nodes_coords])
        self.list_of_decomposed_nodes_coords_Y = np.array([coord[1] for coord in list_of_decomposed_nodes_coords])
        if os.path.exists(self.decomposed_nodes_coords_filename):
            os.remove(self.decomposed_nodes_coords_filename)

        if self.print_hole_geometry_files:
            self.decomposed_nodes_coords_file = open(self.decomposed_nodes_coords_filename, "a")
            for coord in list_of_decomposed_nodes_coords:
                self.decomposed_nodes_coords_file.write(str(coord[1]) + " " + str(-coord[0]) + "\n")
            self.decomposed_nodes_coords_file.close()

    def PenetrationDepthEstimation(self):
        F_th = self.F_th
        V = 4.72e-7 # mm3. Approximate ablated volume for 1 pulses (experimental) and 3W power. For 5 pulses it should be around 2.36e-6
        R_th = self.radius_th
        l_s = V / (np.pi * (0.5 * R_th**2 * (np.log(self.C / F_th)) - 0.25 * self.K * R_th**4))
        return l_s

    # Based on Eq. 27 from Gamaly (2002) - Ablation of solids by femtosecond lasers: Ablation mechanism and ablation thresholds for metals and dielectrics
    '''def EvaporationDepth(self, r): 
        if r >= self.radius_th:
            return 0
        else:
            q = self.C * np.exp(-self.K * r**2)
            d_ev = 0.5 * self.l_s * (np.log(q) - np.log(self.F_th))
            return d_ev'''

    # Based on eq. 6 in Woodfield (2024)
    def EvaporationDepth(self, r):
        if r >= self.r_ast_max:
            return 0.0
        else:
            delta_pen = self.delta_pen
            F_p = self.F_p
            q_ast = self.q_ast
            omega_0 = self.omega_0
            import math
            z_ast = delta_pen * (math.log(F_p / (delta_pen * q_ast)) - 2.0 * (r / omega_0)**2)
            return z_ast

    def CreateResultsFile(self, filename):
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
        with h5py.File(filename, 'a') as f:
            f.attrs['ambient_temperature'] = self.T0
            f.attrs['pulse_energy'] = self.Q
            f.attrs['specific_heat_capacity'] = self.cp
            f.attrs['density'] = self.rho
            f.attrs['conductivity'] = self.conductivity
            # Create a dataset to store the radii
            dataset = f.create_dataset('radii', (self.radii.shape), dtype=self.radii.dtype)
            dataset[:] = self.radii[:]
            f.create_group('temperature_increments')

    def WriteResults(self, filename, process_info):
        step = process_info[KratosMultiphysics.STEP]
        time = step = process_info[KratosMultiphysics.TIME]
        # Open the HDF5 file.
        with h5py.File(filename, 'a') as f:
            assert self.radii.shape == self.temperature_increments.shape
            # Create a dataset to store the radii and temperatures data.
            dataset = f['/temperature_increments'].create_dataset(str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype)
            # Write the radii and temperatures data to the dataset.
            dataset[:] = self.temperature_increments
            # Add a time label to the dataset.
            dataset.attrs["time"] = time

    def _get_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")
        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)
        ## Elements
        ## Note that we check for the elements that require substitution to allow for custom elements
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        element_list = ["EulerianConvDiff","LaplacianElement","MixedLaplacianElement","AdjointHeatDiffusionElement","QSConvectionDiffusionExplicit","DConvectionDiffusionExplicit","LaserAxisymmetricEulerianConvectionDiffusion"]
        if element_name in element_list:
            num_nodes_elements = 0
            if (len(self.main_model_part.Elements) > 0):
                for elem in self.main_model_part.Elements:
                    num_nodes_elements = len(elem.GetNodes())
                    break
            num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
            if not num_nodes_elements:
                num_nodes_elements = domain_size + 1
            name_string = f"{element_name}{domain_size}D{num_nodes_elements}N"
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)
        ## Conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        condition_list = ["FluxCondition","ThermalFace","AxisymmetricThermalFace","LineCondition","SurfaceCondition"]
        if condition_name in condition_list:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break
            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
            if not num_nodes_conditions:
                num_nodes_conditions = domain_size
            name_string = f"{condition_name}{domain_size}D{num_nodes_conditions}N"
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)
        return self.settings["element_replace_settings"]
