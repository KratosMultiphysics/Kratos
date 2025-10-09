import os
import time as timer
from os import environ

from abc import ABC, abstractmethod

import numpy as np


environ["OMP_NUM_THREADS"] = "4" #TODO: why 4? Mkae into a parameter?

import KratosMultiphysics

from KratosMultiphysics.ConvectionDiffusionApplication import (
    convection_diffusion_transient_solver,
)

import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication


if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos


def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolver(model, custom_settings)


class LaserDrillingTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, ABC):
    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the settings in base class
        super().__init__(model, custom_settings)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.delta_time = self.main_model_part.ProcessInfo[
            KratosMultiphysics.DELTA_TIME
        ]  # TODO: Make delta_time a local variable?

        self.jump_between_pulses_counter += self.delta_time
        error_in_delta_time = abs(self.jump_between_pulses_counter - self.time_jump_between_pulses)

        numerical_error = 1e-16  # TODO: Make it global or a parameter?

        # If we just started the simulation or
        # if the time elapsed since the previous pulse (numerically) equals the pulse period
        # we apply a pulse and ablate elements
        if self.pulse_number == 0 or error_in_delta_time < numerical_error:
            self.jump_between_pulses_counter = 0
            self.pulse_number += 1
            self.ResetTemperatureField()

            self.ImposeLaserDeltaTemperature()

            self.RemoveElementsByAblation()
            self.AdjustTemperatureFieldAfterAblation()

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
        spot_diameter = self.beam_waist_diameter * np.sqrt(
            1.0 + ((self.focus_z_offset + self.z_ast_max) / self.rayleigh_length) ** 2
        )
        return spot_diameter

    def ComputePeakFluence(self):
        """
        Computes the peak fluence of a gaussian pulse from its energy and waist radius
        Source: Woodfield 2024, eq (5)

        Parameters
        ----------
        None

        Returns
        -------
        The peak fluence
        """
        return 2.0 * self.Q / (np.pi * self.omega_0**2)  # J/mm2

    def ComputeMaximumAblationRadius(self):
        """
        Computes the maximum ablation radius of a gaussian pulse from its energy and waist radius
        Source: Woodfield 2024, eq (7)

        Parameters
        ----------
        None

        Returns
        -------
        The maximum ablation radius
        """
        return self.omega_0 * np.sqrt(0.5 * np.log(self.F_p / (self.delta_pen * self.q_ast)))

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "compute_vaporisation"             : false,
            "consider_material_refraction"     : false,
            "adjust_T_field_after_ablation"    : false,
            "print_hole_geometry_files"        : false,
            "print_debug_info"                 : false,
            "decomposed_nodes_coords_filename" : "hole_coords_q_ast=q_ast+delta_pen+mesh_type+mesh_size.txt"
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def SetParameters(self):
        # TODO: - Utilitzar GUI de GiD de la LaserDrilling Application per a generar un cas amb 2 materials i
        #         revisar si pots llegir-los utilitzant les funcions del solver base ja existents:
        #         ConvectionDiffusionTransientSolver -> ConvectionDiffusionSolver -> PrepareModelPart() -> self.import_materials()
        #         (applications\ConvectionDiffusionApplication\python_scripts\convection_diffusion_solver.py)
        #         Si no funcinoa, pots provar la GUI de GiD de la ConvectionDiffusion Application per a generar un cas amb 2 materials
        #       - Segurament hauras de generalitzar aquesta funcio per a mes d'1 material: tot el que es llegeixi de self.material_settings
        #         hauria de llegir-se del json de materials (LaserDrillingMaterials.json)

        self.some_elements_are_above_the_evap_temp = False
        self.jump_between_pulses_counter = 0
        self.pulse_number = 0
        self.print_hdf5_and_gnuplot_files = False  # TODO: Make into a parameter

        # Load the material parameters file
        # self.settings is the "solver_settings" section of ProjectParameters.json
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, "r") as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        self.material_settings = materials["properties"][0]["Material"]

        # Load the environment parameters file
        environment_filename = self.settings["material_import_settings"]["environment_filename"].GetString()

        with open(environment_filename, "r") as parameter_file:
            environment = KratosMultiphysics.Parameters(parameter_file.read())

        self.environment_settings = environment["properties"][0]["Material"]

        # Load the laser parameters file
        laser_filename = self.settings["material_import_settings"]["laser_filename"].GetString()

        with open(laser_filename, "r") as parameter_file:
            laser = KratosMultiphysics.Parameters(parameter_file.read())

        self.laser_settings = laser["properties"][0]["Material"]

        """ 
        TODO: instead of hardcoding default values for the parameters, wouldn't it be better if when a parameter is not defined the program failed?
        In this way, we would ensure that the user is sure that they want to set a specific value for a variable. And they
        would not be able to accidentally run a case with values that were not chosen by them
        """
        # Laser parameters
        if not self.laser_settings["Variables"].Has("average_laser_power"):
            self.average_laser_power = 18
        else:
            self.average_laser_power = self.laser_settings["Variables"]["average_laser_power"].GetDouble()

        # pulse frequency = repetition rate
        if not self.laser_settings["Variables"].Has("pulse_frequency"):
            self.pulse_frequency = 2e5
        else:
            self.pulse_frequency = self.laser_settings["Variables"]["pulse_frequency"].GetDouble()

        if not self.laser_settings["Variables"].Has("beam_waist_diameter"):
            self.beam_waist_diameter = 0.0179
        else:
            self.beam_waist_diameter = self.laser_settings["Variables"]["beam_waist_diameter"].GetDouble()

        if not self.laser_settings["Variables"].Has("Rayleigh_length"):
            self.rayleigh_length = 0.409
        else:
            self.rayleigh_length = self.laser_settings["Variables"]["Rayleigh_length"].GetDouble()

        if not self.laser_settings["Variables"].Has("focus_Z_offset"):
            self.focus_z_offset = 0.4
        else:
            self.focus_z_offset = self.laser_settings["Variables"]["focus_Z_offset"].GetDouble()

        if not self.laser_settings["Variables"].Has("reference_T_after_laser"):
            self.reference_T_after_laser = 298.15
        else:
            self.reference_T_after_laser = self.laser_settings["Variables"]["reference_T_after_laser"].GetDouble()

        self.z_ast_max = 0.0  # TODO: Parameter? Remove, since the spot_diameter is unused?

        # self.spot_diameter = self.ComputeSpotDiameter()

        if not self.material_settings["Variables"].Has("VAPORISATION_TEMPERATURE"):
            self.T_e = 1000.0
        else:
            self.T_e = self.material_settings["Variables"]["VAPORISATION_TEMPERATURE"].GetDouble()

        if not self.settings.Has("compute_vaporisation"):
            self.compute_vaporisation = False
        else:
            self.compute_vaporisation = self.settings["compute_vaporisation"].GetBool()

        if not self.material_settings["Variables"].Has("IONIZATION_ALPHA"):
            self.ionization_alpha = 1.0  # 0.95
        else:
            self.ionization_alpha = 1.0  # self.material_settings['Variables']['IONIZATION_ALPHA'].GetDouble()

        if not self.material_settings["Variables"].Has("PENETRATION_DEPTH"):
            self.l_s = 0.002148  # mm.
        else:
            self.l_s = self.material_settings["Variables"]["PENETRATION_DEPTH"].GetDouble()

        if not self.material_settings["Variables"].Has("ABLATION_THRESHOLD"):
            self.F_th = 0.010667  # J/mm2
        else:
            self.F_th = self.material_settings["Variables"]["ABLATION_THRESHOLD"].GetDouble()

        if not self.material_settings["Variables"].Has("THERMAL_DEPTH"):
            self.l_th = 0.0007  # mm.
        else:
            self.l_th = self.material_settings["Variables"]["THERMAL_DEPTH"].GetDouble()

        if not self.material_settings["Variables"].Has("ENTHALPY"):
            self.H_ev = (
                4e5  # J/Kg. Value found on the internet for a given epoxy resin. # TODO: find an actual source ffs
            )
        else:
            self.H_ev = self.material_settings["Variables"]["ENTHALPY"].GetDouble()



        if not self.settings.Has("print_hole_geometry_files"):
            self.print_hole_geometry_files = False
        else:
            self.print_hole_geometry_files = self.settings["print_hole_geometry_files"].GetBool()

        self.Q = self.average_laser_power / self.pulse_frequency  # Energy per pulse
        self.time_jump_between_pulses = 1.0 / self.pulse_frequency  # TODO: rename to something like pulse_period?
        self.cp = self.material_settings["Variables"]["SPECIFIC_HEAT"].GetDouble()
        self.conductivity = self.material_settings["Variables"]["CONDUCTIVITY"].GetDouble()
        self.rho = self.material_settings["Variables"]["DENSITY"].GetDouble()

        # Environment data
        if not self.environment_settings["Variables"].Has("ambient_temperature"):
            self.ambient_temperature = 298.15
        else:
            self.ambient_temperature = self.environment_settings["Variables"]["ambient_temperature"].GetDouble()

        self.T0 = self.ambient_temperature  # The initial temperature is equal to that of the ambient (?)

        self.kappa = self.conductivity / (self.rho * self.cp)
        self.ablation_energy_fraction = self.ionization_alpha
        self.evaporation_energy_fraction = 1.0 - self.ionization_alpha

        ## 2024 Woodfield - Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue
        ## Laser data
        self.omega_0 = 0.5 * self.ComputeSpotDiameter()  # self.R_far # mm # TODO: what is omega_0? omega typically denotes the waist radius but here it is the spot radius?

        y_limit = 2.0 * self.omega_0
        self.hole_theoretical_Y_coords = np.linspace(0.0, float(y_limit), 101)
        self.hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_Y_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)
        self.F_p = self.ComputePeakFluence()

        ## Material calibration using experiments
        if not self.material_settings["Variables"].Has("OPTICAL_PENETRATION_DEPTH"):
            self.delta_pen = 5e-4  # mm
        else:
            self.delta_pen = self.material_settings["Variables"]["OPTICAL_PENETRATION_DEPTH"].GetDouble()

        if not self.material_settings["Variables"].Has("ENERGY_PER_VOLUME_THRESHOLD"):
            self.q_ast = 10.0  # J/mm3
        else:
            self.q_ast = self.material_settings["Variables"]["ENERGY_PER_VOLUME_THRESHOLD"].GetDouble()

        if not self.material_settings["Variables"].Has("REFRACTIVE_INDEX"):
            self.refractive_index_n = 1.5
        else:
            self.refractive_index_n = self.material_settings["Variables"]["REFRACTIVE_INDEX"].GetDouble()

        if not self.settings.Has("consider_material_refraction"):
            self.consider_material_refraction = False
        else:
            self.consider_material_refraction = self.settings["consider_material_refraction"].GetBool()

        if self.material_settings["compute_optical_penetration_depth_using_refractive_index"].GetBool():
            self.l_s = self.ComputeOpticalPenetrationDepth()  
            self.delta_pen = self.l_s

        if self.material_settings["compute_energy_per_unit_volume_threshold_using_enthalpy_and_ionization"].GetBool():
            self.ionizarion_energy_per_volume_threshold = self.ComputeIonizationEnergyPerUnitVolumeThreshold()
            # TODO: typo "ionizaRion"?
            self.use_enthalpy_and_ionization = True
        else:
            self.use_enthalpy_and_ionization = False


        if not self.settings.Has("decomposed_nodes_coords_filename"):
            self.decomposed_nodes_coords_filename = "hole_coords_q_ast=q_ast+delta_pen+mesh_type+mesh_size.txt"
        else:
            self.decomposed_nodes_coords_filename = self.settings["decomposed_nodes_coords_filename"].GetString()

        self.r_ast_max = self.ComputeMaximumAblationRadius()

        if not self.settings.Has("adjust_T_field_after_ablation"):
            self.adjust_T_field_after_ablation = False
        else:
            self.adjust_T_field_after_ablation = self.settings["adjust_T_field_after_ablation"].GetBool()

        if not self.settings.Has("print_debug_info"):
            self.print_debug_info = False
        else:
            self.print_debug_info = self.settings["print_debug_info"].GetBool()

        self.analytical_ablated_volume_in_n_pulses = 0.0

        materials_list = materials["properties"].values() # List with dictionaries that contain the properties of each material
        for material in materials_list:
            full_material_part_name = material["model_part_name"].GetString()
            prefix = "ThermalModelPart."
            material_part_name = full_material_part_name.removeprefix(prefix)

            material_part = self.main_model_part.GetSubModelPart(material_part_name)
            material_settings = material["Material"] # Dictionary with the material's properties
            thermal_energy_per_volume = material_settings["Variables"]["ENERGY_PER_VOLUME_THRESHOLD"].GetDouble()
            for elem in material_part.Elements:
                elem.SetValue(LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME, thermal_energy_per_volume)


        l_th_in_meters = self.l_th * 1e-3
        kappa_in_square_meters = self.kappa * 1e-6
        self.thermal_penetration_time = l_th_in_meters**2 / kappa_in_square_meters


    def ComputeMaximumDepth(self):
        # TODO: Are we sure that what I understand to be the node at the axis of symmetry is
        # always the deepest one? Shouldn't we search for the deepest among the list?
        maximum_depth = self.list_of_ablated_nodes_coords_X[0]
        return maximum_depth

    def ComputeOpticalPenetrationDepth(self):
        # TODO: Find a source for this
        # TODO: Shouldn't this be a parameter read from the JSON parameters?
        light_lambda = 550e-6  # mm, light wavelength

        # TODO: Is there a reason not to combine the following two lines?
        epoxy_n = self.refractive_index_n
        n = epoxy_n
        A = 4.0 * n / ((n + 1) ** 2 + n**2)


        l_s = 0.25 * light_lambda * A / np.pi
        return l_s

    def UpdateLaserRelatedParameters(self):
        self.z_ast_max = self.ComputeMaximumDepth()
        # TODO: Maybe implement a function ComputeWaist that does this? For code clarity and
        # so we don't forget to multiply by 0.5. 
        self.omega_0 = 0.5 * self.ComputeSpotDiameter()
        self.F_p = self.ComputePeakFluence()
        self.r_ast_max = self.ComputeMaximumAblationRadius()


    def AdjustTemperatureFieldAfterAblation(self):
        if not self.adjust_T_field_after_ablation:
            return
        # TODO: check if self.T0 - self.reference_T_after_laser = 0.
        # In that case, new_temperature = old_temperature and this loop can be skipped
        # I don't know if it's worth doing this modification, but right now,
        # self.T0 - self.reference_T_after_laser = 0 (15-04-2025)
        for node in self.main_model_part.Nodes:
            old_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            new_temperature = old_temperature + self.T0 - self.reference_T_after_laser
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temperature)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temperature)
        for elem in self.main_model_part.Elements:
            elem.CalculateOnIntegrationPoints(KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo)

    def ImposeLaserDeltaTemperature(self):
        """
        Calls a function that applies the temperature increase due to the laser pulse.
        Depending on the parameters of the simulation, it decides which function to call.
        For example, it chooses between applying the pulse with or without refraction.
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # TODO: Is there any advantage to writing the conditional in this way, instead of
        # having the clause not be negated?
        if not self.consider_material_refraction:
            self.ImposeTemperatureIncreaseDueToLaser()
        else:
            self.ImposeTemperatureIncreaseDueToLaserWithRefraction()

    def Initialize(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).Initialize()
        self.starting_time = timer.time()

        self.AllocateKratosMemory()

        self.SetParameters()

        # TODO: change python arrays into numpy arrays?
        # TODO: Explain what these are
        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []
        self.list_of_lists_of_decomposed_nodes_X = []
        self.list_of_lists_of_decomposed_nodes_Y = []
        self.full_list_of_ablated_nodes_coords_X = []
        self.full_list_of_ablated_nodes_coords_Y = []



        # TODO: Initial condition, ambient temperature. Do this using GUI!
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)

        self.IdentifyInitialSurfaceNodes()

        self.ResidualHeatStage()




    def IdentifyInitialSurfaceNodes(self):
        # TODO: Why are they called 'decomposed' in this function?
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
            self.main_model_part.CreateSubModelPart("BoundaryPart")
            self.boundary_part = self.main_model_part.GetSubModelPart("BoundaryPart")

        self.boundary_part.AddElements(self.list_of_decomposed_elements_ids)
        self.boundary_part.AddNodes(self.list_of_decomposed_nodes_ids)

    
    def ResetTemperatureField(self):
        # TODO: why would we do this?
        if self.adjust_T_field_after_ablation:
            reference_temp = self.reference_T_after_laser
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, reference_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, reference_temp)

    
    def ComputePulseHoleAndAddToTotalHole(self):
        for i, Y_coord in enumerate(self.hole_theoretical_Y_coords):
            self.hole_theoretical_X_coords[i] += self.EvaporationDepth(Y_coord)

        if self.pulse_number == 1:
            import copy

            self.one_pulse_hole_theoretical_X_coords = copy.deepcopy(self.hole_theoretical_X_coords)
            self.one_pulse_hole_theoretical_Y_coords = self.hole_theoretical_Y_coords

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        
        self.UpdateLaserRelatedParameters()
        self.ComputePulseHoleAndAddToTotalHole()

    def Finalize(self):
        super().Finalize()
        elapsed_time = timer.time() - self.starting_time
        print("\nElapsed_time:", elapsed_time, "\n")
        if self.print_hdf5_and_gnuplot_files:
            self.decomposed_volume_file.close()
            self.temperature_alpha_file.close()
            self.time_alpha_file.close()






   

    def sortSecond(self, val):
        """
        Helper function used, for instance, as a key for sorting a list. Returns the second element of the list.

        Parameters
        ----------
        val: list
             The list.

        Returns
        -------
        The second element of the list
        """

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
        self.list_of_decomposed_nodes_ids = np.array(
            list(set(list_of_decomposed_nodes_ids))
        )  # TODO: why do list(set(x))? To remove duplicates?
        self.list_of_decomposed_elements_ids = np.array(list_of_decomposed_elements_ids)

        if not self.main_model_part.HasSubModelPart("BoundaryPart"):
            self.main_model_part.CreateSubModelPart("BoundaryPart")
            self.boundary_part = self.main_model_part.GetSubModelPart("BoundaryPart")
        else:
            self.main_model_part.RemoveSubModelPart(self.boundary_part)
            self.main_model_part.CreateSubModelPart("BoundaryPart")
            self.boundary_part = self.main_model_part.GetSubModelPart("BoundaryPart")
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

        # TODO: move elsewhere
        # Export a list of the decomposed nodes
        if os.path.exists(self.decomposed_nodes_coords_filename):
            os.remove(self.decomposed_nodes_coords_filename)

        if self.print_hole_geometry_files:
            self.decomposed_nodes_coords_file = open(self.decomposed_nodes_coords_filename, "a")
            for coord in list_of_decomposed_nodes_coords:
                self.decomposed_nodes_coords_file.write(str(coord[1]) + " " + str(-coord[0]) + "\n")
            self.decomposed_nodes_coords_file.close()




    # Based on eq. 6 in Woodfield (2024)
    def EvaporationDepth(self, r):
        """
        Calculates the depth of the ablated cavity as a function of radial position for a single gaussian pulse according to eq. (6) in Woodfield (2024).

        Parameters
        ----------
        r: float
            Radial coordinate with respect to the axis of the beam

        Returns
        -------
        float
            The depth of the cavity at the specified radial coordinate r
        """
        if r >= self.r_ast_max:
            return 0.0
        else:
            delta_pen = self.delta_pen
            F_p = self.F_p
            q_ast = self.q_ast
            omega_0 = self.omega_0
            import math

            z_ast = delta_pen * (math.log(F_p / (delta_pen * q_ast)) - 2.0 * (r / omega_0) ** 2)
            return z_ast


    
    def _get_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2, 3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")
        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)
        ## Elements
        ## Note that we check for the elements that require substitution to allow for custom elements
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        element_list = [
            "EulerianConvDiff",
            "LaplacianElement",
            "MixedLaplacianElement",
            "AdjointHeatDiffusionElement",
            "QSConvectionDiffusionExplicit",
            "DConvectionDiffusionExplicit",
            "LaserAxisymmetricEulerianConvectionDiffusion",
        ]
        if element_name in element_list:
            num_nodes_elements = 0
            if len(self.main_model_part.Elements) > 0:
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
        condition_list = [
            "FluxCondition",
            "ThermalFace",
            "AxisymmetricThermalFace",
            "LineCondition",
            "SurfaceCondition",
        ]
        if condition_name in condition_list:
            num_nodes_conditions = 0
            if len(self.main_model_part.Conditions) > 0:
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break
            num_nodes_conditions = (
                self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
            )
            if not num_nodes_conditions:
                num_nodes_conditions = domain_size
            name_string = f"{condition_name}{domain_size}D{num_nodes_conditions}N"
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)
        return self.settings["element_replace_settings"]

    @abstractmethod
    def ImposeTemperatureIncreaseDueToLaser(self):
        """
        Increases the temperature as an effect of the energy deposition by the laser pulse.
        Does not take into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.
        """
        pass

    @abstractmethod
    def ImposeTemperatureIncreaseDueToLaserWithRefraction(self):
        """
        Increases the temperature as an effect of the energy deposition by the laser pulse.
        Takes into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.
        """
        pass

    @abstractmethod
    def ComputeIonizationEnergyPerUnitVolumeThreshold(self):
        """
        Must be implemented by child classes.
        """
        pass
