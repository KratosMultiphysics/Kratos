import os
import time as timer
from os import environ

from abc import ABC, abstractmethod

import h5py
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import gamma as GammaFunction
from scipy.integrate import quad
from scipy.optimize import root_scalar


# from sympy import *


environ["OMP_NUM_THREADS"] = "4"

import KratosMultiphysics

from KratosMultiphysics.ConvectionDiffusionApplication import (
    convection_diffusion_transient_solver,
)

import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
from KratosMultiphysics.LaserDrillingApplication.fem_tools import SurfaceFEMProjector

if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

from KratosMultiphysics import Logger
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility


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

    #    def ComputeSpotDiameter(self):
    #        """
    #        Computes the 1/e^2 diameter of a gaussian beam at position self.focus_z_offset + self.z_ast_max.
    #        See: https://en.wikipedia.org/wiki/Gaussian_beam#Evolving_beam_width
    #
    #        Parameters
    #        ----------
    #        None
    #
    #        Returns
    #        -------
    #        spot_diameter: float
    #            The diameter of the gaussian beam
    #        """
    #        spot_diameter = self.beam_waist_diameter * np.sqrt(
    #            1.0 + ((self.focus_z_offset + self.z_ast_max) / self.rayleigh_length) ** 2
    #        )
    #        return spot_diameter

    #    def ComputePeakFluenceGaussian(self):
    #        """
    #        Computes the peak fluence of a gaussian pulse from its energy and waist radius
    #        Source: Woodfield 2024, eq (5)
    #
    #        Parameters
    #        ----------
    #        None
    #
    #        Returns
    #        -------
    #        The peak fluence
    #        """
    #        return 2.0 * self.Q / (np.pi * self.omega_0**2)  # J/mm2

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
            self.average_laser_power = 18  # W
        else:
            self.average_laser_power = self.laser_settings["Variables"]["average_laser_power"].GetDouble()

        # pulse frequency = repetition rate
        if not self.laser_settings["Variables"].Has("pulse_frequency"):
            self.pulse_frequency = 2e5  # Pulses per second (Hz)
        else:
            self.pulse_frequency = self.laser_settings["Variables"]["pulse_frequency"].GetDouble()

        self.Q = self.average_laser_power / self.pulse_frequency  # Energy per pulse
        self.time_jump_between_pulses = 1.0 / self.pulse_frequency  # TODO: rename to something like pulse_period?

        if not self.laser_settings["Variables"].Has("beam_waist_diameter"):
            self.beam_waist_diameter = 16.12e-3  # mm
        else:
            self.beam_waist_diameter = self.laser_settings["Variables"]["beam_waist_diameter"].GetDouble()

        ## 2024 Woodfield - Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue
        ## Laser data
        # TODO: this is wrong, check the implications of this mistake. The correct line is the one following
        # self.omega_0 = 0.5 * self.ComputeSpotDiameter()  # self.R_far # mm
        self.omega_0 = self.beam_waist_diameter / 2  # Waist radius (m) of the Gaussian laser spot (Woodfield 2024)

        if not self.laser_settings["Variables"].Has("Rayleigh_length"):
            self.rayleigh_length = 360.4e-3  # mm
        else:
            self.rayleigh_length = self.laser_settings["Variables"]["Rayleigh_length"].GetDouble()

        if not self.laser_settings["Variables"].Has("focus_Z_offset"):
            self.focus_z_offset = 0.4  # mm
        else:
            self.focus_z_offset = self.laser_settings["Variables"]["focus_Z_offset"].GetDouble()

        if not self.laser_settings["Variables"].Has("reference_T_after_laser"):
            self.reference_T_after_laser = 298.15  # K
        else:
            self.reference_T_after_laser = self.laser_settings["Variables"]["reference_T_after_laser"].GetDouble()

        # self.z_ast_max = 0.0  # TODO: Parameter? Remove, since the spot_diameter is unused?

        # Fluence
        table_fluence_types = ["table-not-normalized", "table-normalized"]
        expression_fluence_types = ["expression"]
        self.valid_fluence_types = ["super-gaussian"] + table_fluence_types + expression_fluence_types

        if not self.laser_settings["Variables"].Has("fluence_type"):
            self.fluence_type = "super-gaussian"
            self.gaussian_order = 2
            Logger.PrintWarning("Warning", "Fluence type not specified, defaulting to a super-gaussian of order 2")
        else:
            self.fluence_type = self.laser_settings["Variables"]["fluence_type"].GetString()

        # Choose the fluence function: either a super-gaussian, a function from a table or
        # a function from an expression
        # TODO: documentation: since the simulation is axisymmetric, the pulse in the table has to be
        # symmetric around x=0, peak at 0, and be monotonically decreasing. The table must only contain
        # values for x>=0 (-inf < x < +inf as a 3D coordinate, 0 <= r < +inf when we think of an
        # axisymmetric "slice"). The y values must be nonnegative. The table must be sorted.

        if self.fluence_type == "super-gaussian":
            Logger.PrintInfo("Fluence type", "Using a super-gaussian fluence")
            if not self.laser_settings["Variables"].Has("gaussian_order"):
                self.gaussian_order = 2
                Logger.PrintWarning("Warning", "Gaussian order not specified, defaulting to order 2")
            else:
                self.gaussian_order = self.laser_settings["Variables"]["gaussian_order"].GetInt()

            self.FluenceFunction = self.FluenceSuperGaussian
            self.F_p = self.ComputePeakFluence()

        elif self.fluence_type in table_fluence_types:
            # Check for table filename
            if not self.laser_settings["Variables"].Has("table_filename"):
                Logger.PrintWarning("Warning", "Fluence table filename is not specified")
                raise KeyError("Fluence table filename not specified in laser_settings['Variables']")

            fluence_table_filename = self.laser_settings["Variables"]["table_filename"].GetString()

            # Load and validate table
            try:
                fluence_table_np = np.loadtxt(fluence_table_filename, delimiter=",")
                self.CheckPulseTableProperties(fluence_table_np)
                Logger.PrintInfo("Info", "Fluence table validation succeeded")
            except FileNotFoundError:
                Logger.PrintWarning("Error", f"Fluence table file '{fluence_table_filename}' not found")
                raise
            except (ValueError, TypeError) as e:
                Logger.PrintWarning("Error", f"Fluence table validation failed: {str(e)}")
                raise

            fluence_table_max_r = fluence_table_np[-1, 0]

            self.gaussian_order = None
            self.F_p = None

            if self.fluence_type == "table-not-normalized":
                Logger.PrintInfo(
                    "Fluence type",
                    "Using the fluence function specified as a table (not necessarily with the correct energy)",
                )

                self.FluenceFunction = self.LoadFluenceFromTable(
                    fluence_table_filename, fluence_table_max_r, normalize=True
                )

            elif self.fluence_type == "table-normalized":
                # Assume the table provided by the user is correctly normalized to the total pulse energy
                Logger.PrintInfo(
                    "Fluence type",
                    "Using the fluence function specified as a table with the correct normalization to total pulse energy",
                )

                self.FluenceFunction = self.LoadFluenceFromTable(
                    fluence_table_filename, fluence_table_max_r, normalize=False
                )

            # We are not going to be using the table as a np array
            del fluence_table_np

        elif self.fluence_type in expression_fluence_types:
            Logger.PrintInfo("Fluence type", "Using a fluence function specified as an analytical expression")

            self.gaussian_order = None
            self.F_p = None

            raise NotImplementedError
        else:
            raise ValueError("Invalid fluence_type")

        # Axial energy distribution
        self.axial_energy_distribution_function = self.AxialDistributionBeerLambert

        # self.spot_diameter = self.ComputeSpotDiameter()

        if not self.material_settings["Variables"].Has("VAPORISATION_TEMPERATURE"):
            self.T_e = 1000.0  # K
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
            self.l_s = 0.2148e-3  # mm.
        else:
            self.l_s = self.material_settings["Variables"]["PENETRATION_DEPTH"].GetDouble()

        if not self.material_settings["Variables"].Has("ABLATION_THRESHOLD"):
            self.F_th = 0.010667  # J/mm2 ?
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

        """
        TODO: remove
        if not self.project_parameters["problem_data"].Has("mesh_size"):
            self.mesh_size = "coarse"
        else:
            self.mesh_size = self.project_parameters["problem_data"]["mesh_size"].GetString()

        if not self.project_parameters["problem_data"].Has("mesh_type"):
            self.mesh_type = "unstructured"
        else:
            self.mesh_type = self.project_parameters["problem_data"]["mesh_type"].GetString()
        """

        if not self.settings.Has("print_hole_geometry_files"):
            self.print_hole_geometry_files = False
        else:
            self.print_hole_geometry_files = self.settings["print_hole_geometry_files"].GetBool()

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

        y_limit = 2.0 * self.omega_0
        self.hole_theoretical_Y_coords = np.linspace(0.0, float(y_limit), 101)
        self.hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_Y_coords = np.linspace(0.0, 0.0, 101)
        self.one_pulse_hole_theoretical_X_coords = np.linspace(0.0, 0.0, 101)

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
            self.ComputeOpticalPenetrationDepth()  # TODO: Better to return the value instead of modifying a global?
            self.delta_pen = self.l_s

        if self.material_settings["compute_energy_per_unit_volume_threshold_using_enthalpy_and_ionization"].GetBool():
            self.ionizarion_energy_per_volume_threshold = self.ComputeIonizationEnergyPerUnitVolumeThreshold()
            # TODO: typo "ionizaRion"?
            self.use_enthalpy_and_ionization = True
        else:
            self.use_enthalpy_and_ionization = False

        """
        TODO: see the ToDo list
        self.decomposed_nodes_coords_filename = (
            "hole_coords_q_ast="
            + str(self.q_ast)
            + "_delta_pen="
            + str(self.delta_pen)
            + "_"
            + self.mesh_type
            + "_"
            + self.mesh_size
            + ".txt"
        )
        """

        if not self.settings.Has("decomposed_nodes_coords_filename"):
            self.decomposed_nodes_coords_filename = "hole_coords_q_ast=q_ast+delta_pen+mesh_type+mesh_size.txt"
        else:
            self.decomposed_nodes_coords_filename = self.settings["decomposed_nodes_coords_filename"].GetString()

        if not self.settings.Has("adjust_T_field_after_ablation"):
            self.adjust_T_field_after_ablation = False
        else:
            self.adjust_T_field_after_ablation = self.settings["adjust_T_field_after_ablation"].GetBool()

        if not self.settings.Has("print_debug_info"):
            self.print_debug_info = False
        else:
            self.print_debug_info = self.settings["print_debug_info"].GetBool()

        self.analytical_ablated_volume_in_n_pulses = 0.0

        for properties, i in enumerate(materials["properties"]):
            full_material_part_name_i = i["model_part_name"].GetString()
            prefix = "ThermalModelPart."
            material_part_name_i = full_material_part_name_i.removeprefix(prefix)

            material_part_i = self.main_model_part.GetSubModelPart(material_part_name_i)
            material_settings_i = i["Material"]
            thermal_energy_per_volume_i = material_settings_i["Variables"]["ENERGY_PER_VOLUME_THRESHOLD"].GetDouble()
            for elem in material_part_i.Elements:
                elem.SetValue(LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME, thermal_energy_per_volume_i)

        # Calculated parameters (TODO: separate SetParameters into LoadParameters and CalculateParameters?)
        # Calculate r_ast_max
        self.r_ast_max = self.ComputeMaximumAblationRadius()

        # self.sigma = 0.5 * self.R_far
        # self.K = 1 / (2 * self.sigma**2)
        # import math
        # self.C = self.ablation_energy_fraction * self.Q * self.K / (np.pi * (1 - math.exp(-self.K * self.R_far**2)))
        # self.irradiated_surface_area = np.pi * self.R_far**2

        # if self.ablation_energy_fraction:
        # Find F_th_fraction multiplying F_th so Radius_th = R_far
        # This gives a F_th_fraction of 0.313, a little too flat (maximum not captured)
        # F_th_fraction = self.C * np.exp(-self.K * self.R_far**2) / (self.ablation_energy_fraction * self.Q / A)

        ####################################
        # Calibrated for:                  #
        # residual heat fraction   = 0.05  #
        # Vaporisation temperature = 1000K #
        # Power                    = 3W    #
        ####################################

        # self.F_th = 0.009667 # J/mm2 #F_th_fraction * self.ablation_energy_fraction * self.Q / self.irradiated_surface_area
        # self.radius_th = math.sqrt(math.log(self.C / self.F_th) / self.K)

        # self.l_s = self.PenetrationDepthEstimation() #0.001 # 5 pulses, with evaporation

        # self.l_s = 0.002048
        # self.l_th = 0.33 # micrometers # Thermal depth, assumed

        # l_th_in_meters = self.l_th * 1e-3
        # kappa_in_square_meters = self.kappa * 1e-6
        # self.thermal_penetration_time = l_th_in_meters**2 / kappa_in_square_meters

        # Finite Elements
        # self.n_surface_elements = 10 # number of elements
        # self.sparse_option = True
        # self.surface_nodes_Y_values = np.linspace(0.0, self.R_far, self.n_surface_elements + 1)
        # self.surface_element_size = self.R_far / self.n_surface_elements

        # Debug
        # self.plot_progressive_hole_figures = False

    #    def ComputeMaximumDepth(self):
    #        # TODO: I believe that X[0] is the node at the axis of symmetry.
    #        # Are we sure that X[0] is always the deepest one?
    #        # Shouldn't we search for the deepest among the list?
    #        maximum_depth = self.list_of_ablated_nodes_coords_X[0]
    #        return maximum_depth

    def ComputeOpticalPenetrationDepth(self):
        # TODO: Find a source for this
        # TODO: Shouldn't this be a parameter read from the JSON parameters?
        light_lambda = 550e-6  # mm, light wavelength

        # TODO: Is there a reason not to combine the following two lines?
        epoxy_n = self.refractive_index_n
        n = epoxy_n
        A = 4.0 * n / ((n + 1) ** 2 + n**2)

        # TODO: Why modify the global self.l_s and also return its value?
        self.l_s = 0.25 * light_lambda * A / np.pi
        return self.l_s

    def ComputePeakFluenceSuperGaussian(self):
        """
        Computes the peak fluence of a super-gaussian pulse from its energy and waist radius.
        The gaussian pulse is a super-gaussian pulse of order 2.

        A super-gaussian pulse has intensity I(r) = I_peak * exp(-2 (r/w_0)^n).  The total
        energy in the pulse is the integral of I(r) from 0 to +inf, which can be evaluated
        using https://en.wikipedia.org/wiki/Gaussian_integral#Relation_to_the_gamma_function.

        Parameters
        ----------
        None

        Returns
        -------
        The peak fluence
        """
        b = self.gaussian_order
        F_p = self.Q * b * 2 ** (2 / b) / (2 * np.pi * self.omega_0**2 * GammaFunction(2 / b))
        return F_p  # J/mm2

    def ComputePeakFluence(self):
        if self.fluence_type != "super-gaussian":
            # Think about this because not all distributions are peaked.
            self.F_p = None
        else:
            # By default, use a super-gaussian pulse
            peak_fluence = self.ComputePeakFluenceSuperGaussian()

        return peak_fluence

    def ComputeSuperGaussianMaximumAblationRadius(self, parameters=None):
        """
        Computes the maximum radius of the ablated cavity after one pulse,
        which corresponds to the radius at the surface of the sample.
        See: Woodfield 2024 eq. (7),

        Parameters
        ----------
        None

        Returns
        -------
        r_ast_max: float
            The radius of the cavity at the surface
        """
        n = self.gaussian_order
        r_ast_max = (
            self.omega_0 / np.power(2, 1.0 / n) * np.power(np.log(self.F_p / (self.delta_pen * self.q_ast)), 1.0 / n)
        )
        return r_ast_max

    def FindPulseRoot(self, F, y, parameters=None, numerical_zero=1e-5):
        """
        Finds the root such that F(root) = y for a pulse function F(x), which is symmetric
        around x=mu, peaks at mu, and is monotonically decreasing for x > mu. The function F
        may take no parameters (F(x)) or a parameters dictionary (F(x, parameters)).

        Parameters
        -----------
        F : callable
            The pulse function, either F(x) or F(x, parameters), taking a numeric input x.
        y : float
            The target value such that F(root) = y.
        parameters : dict or None, optional
            Dictionary of parameters for F (e.g., {"mu": 0, "sigma": 1}). If None, F is called
            without parameters (default: None).
        numerical_zero : float, optional
            Threshold for checking if F(root) is close to y (default: 1e-15).

        Returns
        --------
        float
            The value root such that F(root) = y, with root >= mu.

        Raises:
        -------
        ValueError
            If F is not callable, parameters is invalid, y is invalid, no suitable upper bound
            is found, or root-finding fails.
        """
        if not callable(F):
            raise ValueError("F must be a callable function")

        if parameters is not None and not isinstance(parameters, dict):
            raise ValueError("parameters must be a dictionary or None")

        if y <= 0:
            raise ValueError("y must be positive since F(x) is assumed positive")

        # Get mu (default to 0 if parameters is None or mu not provided)
        mu = 0 if parameters is None else parameters.get("mu", 0)

        # Evaluate F at peak (x=mu) to test signature and check y
        try:
            F_peak = F(mu) if parameters is None else F(mu, parameters)
        except TypeError as e:
            raise ValueError(f"F evaluation at mu failed: {str(e)}")

        if F_peak < y:
            raise ValueError(f"y={y} is greater than F at the peak, F({mu=})={F_peak}, no solution exists")

        # Find x_max such that F(x_max) < y
        x_max = mu + 0.2  # Initial guess
        max_iterations = 100
        step_factor = 2.0  # Increase x_max exponentially

        for _ in range(max_iterations):
            try:
                F_x = F(x_max) if parameters is None else F(x_max, parameters)
            except TypeError as e:
                raise ValueError(f"F evaluation at x={x_max} failed: {str(e)}")

            if F_x < y:
                break
            if not np.isfinite(F_x):
                raise ValueError(f"F(x={x_max}) is {F_x}")
            x_max = mu + (x_max - mu) * step_factor
        else:
            raise ValueError("Could not find x where F(x) < y after maximum iterations")

        # Check if F(x_max) is valid
        if not np.isfinite(F_x):
            raise ValueError(f"F(x={x_max}) is {F_x}")

        # Define the objective function based on parameters
        if parameters is None:

            def Objective(x):
                return F(x) - y
        else:

            def Objective(x):
                return F(x, parameters) - y

        # Solve using root_scalar with Brent's method
        try:
            result = root_scalar(Objective, bracket=[mu, x_max], method="brentq")
            if not result.converged:
                raise ValueError("Root-finding did not converge")
            root = result.root

            # Verify solution
            F_root = F(root) if parameters is None else F(root, parameters)
            if abs(F_root - y) > numerical_zero:
                raise ValueError("Solution does not satisfy F(root) = y within tolerance")

            return root
        except ValueError as e:
            raise ValueError(f"Failed to find root: {str(e)}")

    def ComputeMaximumAblationRadius(self, parameters=None):
        """
        Computes the maximum ablation radius for a material under a pulsed fluence.

        The ablation radius is the maximum radial distance where the fluence is enough to cause ablation.
        The fluence is modeled based on `fluence_type`, which can be a super-Gaussian (analytical solution),
        a table-based function, or a mathematical expression (numerical root-finding via `FindPulseRoot`).

        Parameters
        ----------
        parameters : dict, optional
            Dictionary of parameters for the fluence function (e.g., {"F_p": 1.0, "gamma": 1.0, "mu": 0}).
            Required if `fluence_type` is "table" or "expression" and the fluence function expects parameters,
            or if `ComputeSuperGaussianMaximumAblationRadius` requires parameters.
            Default is None, used for parameterless fluence functions.

        Returns
        -------
        r_ast_max : float
            The maximum ablation radius.

        Raises
        ------
        ValueError
            If `fluence_type` is invalid, `delta_pen` or `q_ast` are not positive, `FluenceFunction` is not callable,
            the computed radius is negative or non-finite, or `y` is invalid.
        TypeError
            If `FluenceFunction` or `parameters` are incorrectly specified when calling `FindPulseRoot`
            or `ComputeSuperGaussianMaximumAblationRadius`.

        Notes
        -----
        - For `fluence_type="super-gaussian"`, the radius is computed analytically using
        `ComputeSuperGaussianMaximumAblationRadius`.
        - For `fluence_type="table"` or `"expression"`, the radius is found numerically by solving
        `FluenceFunction(r, parameters) = delta_pen * q_ast` using `FindPulseRoot`.
        - The fluence function is expected to be a pulse function (symmetric, unimodal, peaking at
        `mu`, and monotonically decreasing for `r > mu`), consistent with Woodfield (2024).
        - `delta_pen` and `q_ast` are instance attributes defining the penetration coefficient and
        ablation threshold energy, respectively.

        See Also
        --------
        FindPulseRoot : Numerical root-finding for pulse functions.
        ComputeSuperGaussianMaximumAblationRadius : Analytical solution for super-Gaussian fluence.
        """
        # Validate fluence_type
        if not isinstance(self.fluence_type, str):
            raise ValueError("fluence_type must be a string")
        if self.fluence_type not in self.valid_fluence_types:
            raise ValueError(f"fluence_type must be one of {self.valid_fluence_types}")

        # Validate instance attributes
        if not hasattr(self, "delta_pen") or not isinstance(self.delta_pen, (int, float)) or self.delta_pen <= 0:
            raise ValueError("delta_pen must be a positive number")
        if not hasattr(self, "q_ast") or not isinstance(self.q_ast, (int, float)) or self.q_ast <= 0:
            raise ValueError("q_ast must be a positive number")

        # Compute ablation limit
        y = self.delta_pen * self.q_ast
        if not np.isfinite(y) or y <= 0:
            raise ValueError(f"Computed ablation limit y={y} is invalid")

        # Compute ablation radius
        if self.fluence_type == "super-gaussian":
            try:
                r_ast_max = self.ComputeSuperGaussianMaximumAblationRadius(parameters=parameters)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Failed to compute super-Gaussian ablation radius: {str(e)}")
        else:  # table or expression
            if not hasattr(self, "FluenceFunction") or not callable(self.FluenceFunction):
                raise ValueError("FluenceFunction must be a callable representing a table or an expression")
            try:
                # We wrap self.FluenceFunction to adapt it to the signature of F in FindPulseRoot
                def FluenceFunctionWrapped(x, parameters=None):
                    position = {"r": x}
                    return self.FluenceFunction(position, parameters)

                r_ast_max = self.FindPulseRoot(FluenceFunctionWrapped, y, parameters=parameters)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Failed to find root for FluenceFunction: {str(e)}")

        # Validate computed radius
        if not np.isfinite(r_ast_max):
            raise ValueError(f"Computed ablation radius r_ast_max={r_ast_max} is non-finite")
        if r_ast_max < 0:
            raise ValueError(f"Computed ablation radius r_ast_max={r_ast_max} is negative")

        return r_ast_max

    def CheckPulseTableProperties(self, arr):
        """
        Validates that a NumPy array representing (x, y) pairs satisfies specific properties for a pulse
        function table: all x values (first column) are non-negative, all y values (second column) are
        non-negative, and y values are monotonically decreasing.

        Parameters
        ----------
        arr : ndarray
            2D NumPy array of shape (n, m) where n >= 1 (number of pairs) and m >= 2 (at least x and y
            columns). The first column contains x values (e.g., radial distances), and the second column
            contains y values (e.g., fluence or energy).

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If arr is not a 2D NumPy array, has fewer than 2 columns, is empty, contains non-finite values,
            or violates any of the required properties (negative x, negative y, non-strictly increasing x,
            or non-monotonically decreasing y).
        TypeError
            If arr is not a NumPy array.

        Notes
        -----
        - Assumes arr represents a pulse function table, such as fluence vs. radial distance
        - x values must be strictly increasing (x[i] < x[i+1]) to ensure the table is sorted.
        - y values must be non-strictly monotonically decreasing (y[i] >= y[i+1]).

        """

        # Validate input type
        if not isinstance(arr, np.ndarray):
            raise TypeError("arr must be a NumPy array")

        # Check array dimensions
        if arr.ndim != 2:
            raise ValueError("arr must be a 2D array")
        if arr.shape[0] < 1:
            raise ValueError("arr must have at least one row")
        if arr.shape[1] < 2:
            raise ValueError("arr must have at least two columns (x and y)")

        # Check for finite values
        if not np.all(np.isfinite(arr)):
            raise ValueError("arr contains non-finite values (NaN or Inf)")

        # Extract x and y
        x = arr[:, 0]
        y = arr[:, 1]

        # Check x >= 0
        if not np.all(x >= 0):
            raise ValueError("All x values must be non-negative")

        # Check x is strictly increasing (x[i] < x[i+1])
        if arr.shape[0] > 1:  # Only check if more than one row
            if not np.all(np.diff(x) > 0):
                raise ValueError("x values must be strictly increasing (x[i] < x[i+1])")

        # Check y >= 0
        if not np.all(y >= 0):
            raise ValueError("All y values must be non-negative")

        # Check y is monotonically decreasing (y[i] >= y[i+1])
        if arr.shape[0] > 1:  # Only check if more than one row
            if not np.all(np.diff(y) <= 0):
                raise ValueError("y values must be monotonically decreasing (y[i] >= y[i+1])")

    def UpdateLaserRelatedParameters(self):
        # self.z_ast_max = self.ComputeMaximumDepth()
        # self.omega_0 = 0.5 * self.ComputeSpotDiameter()
        # self.F_p = self.ComputePeakFluence()
        # self.r_ast_max = self.ComputeMaximumAblationRadius()
        pass

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
        self.results_filename = "results.h5"
        self.CreateResultsFile(self.results_filename)
        self.temperature_increments = np.array(
            [node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes]
        )

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

        # TODO: change python arrays into numpy arrays? Into a modelpart?
        # TODO: Explain what these are
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
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)

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
        # TODO: Why are the nodes called 'decomposed' in this function?
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

    def StorePreEvaporationTemperature(self):
        for node in self.main_model_part.Nodes:
            pre_evaporation_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            node.SetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE, pre_evaporation_temperature)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)
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

    #    def ResidualHeatStage(self):
    #        """
    #        TODO: Currently, unused. It is overriden by LaserDrillingTransienSolverAblationPlusThermal.ResidualHeatStage
    #        Idea of what it does (WIP): after having applied the laser pulse and having ablated those elements whose energy
    #        is higher than a threshold (Woodfield model), the remaining elements are hot and need to thermalise. Woodfield simply
    #        propagates the heat, regardless of whether the temperature of an element is above or below the vaporisation temperature.
    #        In contrast, this function finds elements where the temperature is above the vaporisation temperature of the material
    #        and evaporates them.
    #        """
    #        if self.evaporation_energy_fraction:
    #            self.projector = SurfaceFEMProjector(
    #                self.n_surface_elements, self.R_far, self.sparse_option
    #            )  # , delta_coefficients)
    #            self.q_interp = self.projector.InterpolateFunctionAndNormalize(self.EnergyPerUnitArea1D)  # , 1.0)
    #            if self.compute_vaporisation:
    #                self.first_evaporation_stage_done = False
    #                self.max_vaporisation_layers = 50
    #                # if self.pulse_number == 1:
    #                #     self.max_vaporisation_layers = 1
    #                self.vaporisation_layer_number = 1
    #                self.some_elements_are_above_the_evap_temp = True
    #                print("Removing elements by evaporation...")
    #                print("Pulse number:", self.pulse_number)
    #                self.last_evaporation_layer_applied = False
    #                self.StorePreEvaporationTemperature()
    #                self.ablation_only_legend_added = False
    #
    #                while self.some_elements_are_above_the_evap_temp:
    #                    if self.vaporisation_layer_number > self.max_vaporisation_layers:
    #                        if not self.max_vaporisation_layers:
    #                            self.ImposeTemperatureIncreaseDueTo1DConduction()
    #                        print("******************************MAXIMUM ITERATIONS EXCEEDED!!!")
    #                        break
    #                    self.ImposeTemperatureIncreaseDueTo1DConduction()
    #                    self.RemoveElementsByEvaporation()
    #
    #                    print("Done!")
    #                else:
    #                    self.ImposeTemperatureIncreaseDueTo1DConduction()

    #    def RemoveElementsByEvaporation(self):
    #        evap_elements_centers_Y = []
    #        evap_elements_volumes = []
    #        delta_temp_elements = []
    #        uncapped_delta_temp_elements = []
    #        self.some_elements_are_above_the_evap_temp = False
    #        for elem in self.boundary_part.Elements:
    #            if elem.Is(KratosMultiphysics.ACTIVE):
    #                temp = elem.CalculateOnIntegrationPoints(
    #                    KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo
    #                )
    #                element_temperature = temp[0]
    #                if element_temperature > self.T_e:
    #                    pre_evap_temp = elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE)
    #                    delta_temp = self.T_e - pre_evap_temp
    #                    uncapped_delta_temp = element_temperature - pre_evap_temp
    #                    delta_temp_elements.append(delta_temp)
    #                    uncapped_delta_temp_elements.append(uncapped_delta_temp)
    #                    self.some_elements_are_above_the_evap_temp = True
    #                    elem.Set(KratosMultiphysics.ACTIVE, False)
    #                    vol = elem.CalculateOnIntegrationPoints(
    #                        LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo
    #                    )
    #                    elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, vol[0])
    #                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
    #                    Y_centroid = elem.GetGeometry().Center().Y
    #                    evap_elements_centers_Y.append(Y_centroid)
    #                    evap_elements_volumes.append(element_volume)
    #                    for node in elem.GetNodes():
    #                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)
    #
    #        number_of_problematic_elements = 0
    #        for elem in self.main_model_part.Elements:
    #            if elem.Is(KratosMultiphysics.ACTIVE):
    #                number_of_decomposed_nodes = 0
    #                temp = elem.CalculateOnIntegrationPoints(
    #                    KratosMultiphysics.TEMPERATURE, self.main_model_part.ProcessInfo
    #                )
    #                element_temperature = temp[0]
    #                for node in elem.GetNodes():
    #                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
    #                        number_of_decomposed_nodes += 1
    #                if number_of_decomposed_nodes == 3:
    #                    elem.Set(KratosMultiphysics.ACTIVE, False)
    #                    elem.CalculateOnIntegrationPoints(
    #                        LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo
    #                    )
    #                    element_volume = elem.GetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME)
    #                    Y_centroid = elem.GetGeometry().Center().Y
    #                    pre_evap_temp = elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE)
    #                    delta_temp = self.T_e - pre_evap_temp
    #                    # delta_temp_elements.append(delta_temp)
    #                    uncapped_delta_temp = element_temperature - pre_evap_temp
    #                    # uncapped_delta_temp_elements.append(uncapped_delta_temp)
    #                    # evap_elements_centers_Y.append(Y_centroid)
    #                    # evap_elements_volumes.append(element_volume)
    #                    number_of_problematic_elements += 1
    #        print("Number_of_problematic_elements:", number_of_problematic_elements)
    #        evap_elements_centers_Y = np.array(evap_elements_centers_Y)
    #        evap_elements_volumes = np.array(evap_elements_volumes)
    #        self.evap_elements_centers_Y = evap_elements_centers_Y[evap_elements_centers_Y.argsort()]
    #        self.evap_elements_volumes = evap_elements_volumes[evap_elements_centers_Y.argsort()]
    #        # print("self.evap_elements_centers_Y:", self.evap_elements_centers_Y)
    #        # print("self.evap_elements_volumes:", self.evap_elements_volumes)
    #
    #        delta_temp_elements = np.array(delta_temp_elements)
    #        self.delta_temp_elements = delta_temp_elements[evap_elements_centers_Y.argsort()]
    #        uncapped_delta_temp_elements = np.array(uncapped_delta_temp_elements)
    #        self.uncapped_delta_temp_elements = uncapped_delta_temp_elements[evap_elements_centers_Y.argsort()]
    #
    #        # Total enthalpy: Energy consumed for the material to vaporize plus the energy for heating up the material to the vaporization temperature.
    #        # Equation (8) in Wang, 2019. 'Thermal effect of femtosecond laser polystyrene processing'
    #        self.evap_elements_enthalpies = (
    #            self.evap_elements_volumes * self.rho * (self.H_ev + self.cp * self.delta_temp_elements)
    #        )
    #
    #        # print("evap_elements_volumes:", evap_elements_volumes)
    #        # print("delta_temp_elements:", delta_temp_elements)
    #        # print("self.evap_elements_enthalpies:", self.evap_elements_enthalpies)
    #
    #        self.support_elements = [[] for i in range(self.n_surface_elements + 1)]
    #
    #        if not self.some_elements_are_above_the_evap_temp:
    #            return
    #
    #        import matplotlib.pyplot as plt
    #
    #        figure, axis = plt.subplots(2, 2, figsize=(15, 8))
    #        label_size = 13
    #        numbers_size = 11
    #        axis[0][0].grid()
    #        axis[0][0].plot(self.evap_elements_centers_Y, self.evap_elements_volumes, color="red", marker="+")
    #        axis[0][0].set_xlabel("radius (mm)")
    #        axis[0][0].set_ylabel("Volumes (mm3)", fontsize=label_size)
    #        axis[1][0].grid()
    #        axis[1][0].plot(self.evap_elements_centers_Y, self.delta_temp_elements, color="blue", marker="o")
    #        # axis[1][0].set_ylim(bottom=0, top=None)
    #        axis[1][0].set_xlabel("radius (mm)")
    #        axis[1][0].set_ylabel("Delta temps (K)", fontsize=label_size)
    #        axis[0][1].grid()
    #        axis[0][1].plot(self.evap_elements_centers_Y, self.evap_elements_enthalpies, color="black", marker="x")
    #        axis[0][1].set_xlabel("radius (mm)")
    #        axis[0][1].set_ylabel("Enthalpies (J)", fontsize=label_size)
    #        axis[1][1].grid()
    #        axis[1][1].plot(self.evap_elements_centers_Y, self.uncapped_delta_temp_elements, color="green", marker="*")
    #        # axis[1][1].set_ylim(bottom=0, top=None)
    #        axis[1][1].set_xlabel("radius (mm)")
    #        axis[1][1].set_ylabel("Uncapped delta temps (K)", fontsize=label_size)
    #        # figure.show()
    #
    #        print("\nEvaporating layer number", self.vaporisation_layer_number, "...")
    #        self.vaporisation_layer_number += 1
    #
    #        if not self.sparse_option:
    #            self.projector.FillUpMassMatrix()
    #        else:
    #            self.projector.FillUpSparseMassMatrix()
    #        self.projector.AssignDeltasToTestFunctionSupports(self.evap_elements_centers_Y, self.support_elements)
    #        self.projector.FillUpDeltasRHS(
    #            self.evap_elements_centers_Y, self.support_elements, self.evap_elements_enthalpies
    #        )
    #        self.u = self.projector.Project()
    #
    #        total_energy = self.projector.CalculateEnergyOfFEMFunction(self.u)
    #
    #        print("Total energy expected =", sum(self.evap_elements_enthalpies))
    #        print("Total energy calculated =", total_energy)
    #
    #        import matplotlib.pyplot as plt
    #
    #        _, axis = plt.subplots(1, 2, figsize=(20, 5))
    #        axis[0].grid()
    #        axis[0].plot(self.projector.X, self.q_interp, color="red", marker="+")
    #        axis[0].plot(self.projector.X, self.u, color="blue", marker="o")
    #
    #        self.q_interp -= self.u
    #
    #        # TODO: rethink this!
    #        for i, q in enumerate(self.q_interp):
    #            if q < 0.0:
    #                self.q_interp[i] = 0.0
    #
    #        axis[0].tick_params(axis="both", which="major", labelsize=numbers_size)
    #        axis[0].tick_params(axis="both", which="minor", labelsize=numbers_size)
    #        plt.subplots_adjust(left=0.1, bottom=None, right=0.9, top=None, wspace=0.22, hspace=None)
    #        axis[0].plot(self.projector.X, self.q_interp, color="black", marker="*")
    #        axis[0].set_ylim(bottom=0.0, top=0.0036)
    #        axis[0].legend(
    #            ["fluence (interpolated)", "fluence (lost)", "fluence (remaining)"],
    #            loc="upper right",
    #            fontsize=label_size,
    #        )
    #        axis[0].set_xlabel("radius (mm)", fontsize=label_size)
    #        axis[0].set_ylabel("Energies (J/mm2)", fontsize=label_size)
    #
    #        self.AddDecomposedNodesToSurfaceList()
    #        self.first_evaporation_stage_done = True
    #
    #        X = self.list_of_decomposed_nodes_coords_X
    #        Y = self.list_of_decomposed_nodes_coords_Y
    #
    #        if not self.ablation_only_legend_added:
    #            if self.ablation_energy_fraction:
    #                self.full_list_of_ablated_nodes_coords_X.append(self.list_of_ablated_nodes_coords_X)
    #                self.full_list_of_ablated_nodes_coords_Y.append(self.list_of_ablated_nodes_coords_Y)
    #            else:
    #                self.full_list_of_ablated_nodes_coords_X.append(0.0 * Y)
    #                self.full_list_of_ablated_nodes_coords_Y.append(Y)
    #            self.ablation_only_legend_added = True
    #
    #        self.list_of_lists_of_decomposed_nodes_X.append(X)
    #        self.list_of_lists_of_decomposed_nodes_Y.append(Y)
    #
    #        axis[1].grid()
    #        list_of_legends = []
    #        p = 1
    #        axis[1].tick_params(axis="both", which="major", labelsize=numbers_size)
    #        axis[1].tick_params(axis="both", which="minor", labelsize=numbers_size)
    #        for list_X, list_Y in zip(self.full_list_of_ablated_nodes_coords_X, self.full_list_of_ablated_nodes_coords_Y):
    #            axis[1].plot(list_Y, -list_X, color="black")
    #            axis[1].set_ylim(bottom=-0.0023, top=None)
    #            converted_num = "Ablation number #" + str(p)
    #            list_of_legends.append(converted_num)
    #            p += 1
    #        i = 1
    #        for list_X, list_Y in zip(self.list_of_lists_of_decomposed_nodes_X, self.list_of_lists_of_decomposed_nodes_Y):
    #            axis[1].plot(list_Y, -list_X)
    #            axis[1].set_ylim(bottom=-0.0023, top=0)
    #            converted_num = "Evap. layer #" + str(i)
    #            list_of_legends.append(converted_num)
    #            i += 1
    #        axis[1].legend(list_of_legends, loc="lower right", fontsize=11)
    #        axis[1].set_xlabel("radius (mm)", fontsize=label_size)
    #        axis[1].set_ylabel("Ablation plus evaporation (mm)", fontsize=label_size)
    #
    #        if self.plot_progressive_hole_figures:
    #            plt.show()
    #
    #        if self.vaporisation_layer_number <= self.max_vaporisation_layers:
    #            self.RetrievePreEvaporationTemperatureState()

    #    def RetrievePreEvaporationTemperatureState(self):
    #        for node in self.main_model_part.Nodes:
    #            node.SetSolutionStepValue(
    #                KratosMultiphysics.TEMPERATURE,
    #                node.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE),
    #            )
    #            node.SetSolutionStepValue(
    #                KratosMultiphysics.TEMPERATURE,
    #                1,
    #                node.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE),
    #            )
    #        for elem in self.main_model_part.Elements:
    #            elem.SetValue(
    #                KratosMultiphysics.TEMPERATURE,
    #                elem.GetValue(LaserDrillingApplication.PRE_EVAPORATION_TEMPERATURE),
    #            )

    #    def ImposeTemperatureIncreaseDueTo1DConduction(self):
    #        """
    #        TODO: I think this is currently unused
    #        """
    #        X = self.list_of_decomposed_nodes_coords_X
    #        Y = self.list_of_decomposed_nodes_coords_Y
    #        # TODO: why these values?
    #        self.minimum_characteristic_Z = 1e6
    #        self.maximum_characteristic_Z = -1e6
    #
    #        for node in self.main_model_part.Nodes:
    #            radius = node.Y
    #            """ if not self.ablation_energy_fraction:
    #                distance_to_surface = 0.0
    #            else: """
    #            F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
    #            distance_to_surface = F(radius)
    #            z = node.X - distance_to_surface
    #            if radius <= self.r_ast_max:
    #                delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
    #                old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
    #                new_temp = old_temp + delta_temp
    #                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
    #                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)
    #        # print("\nResidual heat fraction:", self.ionization_alpha)
    #        # print("\nMaximum characteristic depth:", self.maximum_characteristic_Z)
    #        # problem_characteristic_time_minimum_depth = self.minimum_characteristic_Z**2 / self.kappa
    #        # problem_characteristic_time_maximum_depth = self.maximum_characteristic_Z**2 / self.kappa
    #        # minimum_time_step_for_minimum_depth = 0.1 * problem_characteristic_time_minimum_depth
    #        # minimum_time_step_for_maximum_depth = 0.1 * problem_characteristic_time_maximum_depth
    #        # print("\nThermal problem characteristic time for maximum depth:", problem_characteristic_time_maximum_depth)
    #        # print("\nNecessary time step for maximum depth:", minimum_time_step_for_maximum_depth, '\n')

    #    def EnergyPerUnitArea1D(self, radius):
    #        C = (1 - self.ablation_energy_fraction) * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
    #        q = C * np.exp(-self.K * radius**2)
    #        return q

    #    def InitialThermalConductionTime(self, radius):
    #        # TODO: This is never called. Also, it has radius as a parameter but it is unused. In addition, it does nothing, it just returns a variable that is calculated in SetParameters
    #
    #        # This function returns the characteristic time required for the initial heat distribution from the surface to the interior
    #        # 4.0 (2**2) in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
    #        # C = 4.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * (self.T_e - self.T0)**2)
    #        # t = C * self.EnergyPerUnitArea1D(radius)**2
    #        return self.thermal_penetration_time

    #    def TemperatureVariationInZDueToLaser1D(self, radius, z):
    #        # q = self.EnergyPerUnitArea1D(radius)
    #        q = self.projector.EvaluateFEMFunction(self.q_interp, radius)
    #        t_penetration = self.InitialThermalConductionTime(radius)
    #        # 2.0 in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
    #        C = 2.0 * q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_penetration))
    #        delta_temp = C * np.exp(-(z**2) / (4.0 * self.kappa * t_penetration))
    #        characteristic_Z = np.sqrt(4.0 * self.kappa * t_penetration)
    #        if characteristic_Z <= self.minimum_characteristic_Z:
    #            self.minimum_characteristic_Z = characteristic_Z
    #        if characteristic_Z >= self.maximum_characteristic_Z:
    #            self.maximum_characteristic_Z = characteristic_Z
    #        return delta_temp

    def SolveSolutionStep(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).SolveSolutionStep()
        # TODO: Perhaps the following fits better in FinalizeSolutionStep?
        if self.print_hdf5_and_gnuplot_files:
            self.temperature_increments = np.array(
                [node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes]
            )
            self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    #    def ComputePulseHoleAndAddToTotalHole(self):
    #        for i, Y_coord in enumerate(self.hole_theoretical_Y_coords):
    #            self.hole_theoretical_X_coords[i] += self.EvaporationDepth(Y_coord)
    #
    #        if self.pulse_number == 1:
    #            import copy
    #
    #            self.one_pulse_hole_theoretical_X_coords = copy.deepcopy(self.hole_theoretical_X_coords)
    #            self.one_pulse_hole_theoretical_Y_coords = self.hole_theoretical_Y_coords

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.print_hdf5_and_gnuplot_files:
            decomp_vol = self.MonitorDecomposedVolume()
            decomp_vol *= 1e9  # To convert mm3 into um3
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
        # self.ComputePulseHoleAndAddToTotalHole()

    def Finalize(self):
        super().Finalize()
        elapsed_time = timer.time() - self.starting_time
        print("\nElapsed_time:", elapsed_time, "\n")
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
        plt.plot(x, y)
        plt.show()

    def MonitorDecomposedVolume(self):
        """
        Tallies up the volume of the decomposed elements

        Parameters
        ----------
        None

        Returns
        -------
        The sum of DECOMPOSED_ELEMENTAL_VOLUME
        """
        decomposed_volume = 0.0
        for elem in self.main_model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(
                LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, self.main_model_part.ProcessInfo
            )
            decomposed_volume += out[0]
        return decomposed_volume

    def MonitorEnergy(self):
        """
        Tallies up the energy of the (not decomposed) elements

        Parameters
        ----------
        None

        Returns
        -------
        The sum of THERMAL_ENERGY over all ACTIVE elements
        """
        energy = 0.0
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
                out = elem.CalculateOnIntegrationPoints(
                    LaserDrillingApplication.THERMAL_ENERGY, self.main_model_part.ProcessInfo
                )
                energy += out[0]
        return energy

    def RemoveElementsByAblation(self):
        # TODO: make this method abstract? It is overridden by LaserDrillingTransientSolverAblationPlusThermal.RemoveElementsByAblation
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
                if DeltaX <= d_ev:  # and Y_centroid <= self.radius_th:
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
            """print('\nR_far:', self.R_far)
            print('\nRadius_th:', self.radius_th)
            print("\nDecomposed volume:", self.MonitorDecomposedVolume())"""

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

    def PenetrationDepthEstimation(self):
        # TODO: this function is never called. Should we remove it?
        # TODO: where does this formula come from?
        F_th = self.F_th
        # TODO: make this into a parameter
        V = 4.72e-7  # mm3. Approximate ablated volume for 1 pulses (experimental) and 3W power. For 5 pulses it should be around 2.36e-6
        R_th = self.radius_th
        l_s = V / (np.pi * (0.5 * R_th**2 * np.log(self.C / F_th) - 0.25 * self.K * R_th**4))
        return l_s

    # Based on Eq. 27 from Gamaly (2002) - Ablation of solids by femtosecond lasers: Ablation mechanism and ablation thresholds for metals and dielectrics
    # TODO: remove?
    """def EvaporationDepth(self, r): 
        if r >= self.radius_th:
            return 0
        else:
            q = self.C * np.exp(-self.K * r**2)
            d_ev = 0.5 * self.l_s * (np.log(q) - np.log(self.F_th))
            return d_ev"""

    # Based on eq. 6 in Woodfield (2024)
    # TODO: rename to AblationDepth?
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

            z_ast = delta_pen * (np.log(F_p / (delta_pen * q_ast)) - 2.0 * (r / omega_0) ** 2)
            return z_ast

    def CreateResultsFile(self, filename):
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
        with h5py.File(filename, "a") as f:
            f.attrs["ambient_temperature"] = self.T0
            f.attrs["pulse_energy"] = self.Q
            f.attrs["specific_heat_capacity"] = self.cp
            f.attrs["density"] = self.rho
            f.attrs["conductivity"] = self.conductivity
            # Create a dataset to store the radii
            dataset = f.create_dataset("radii", (self.radii.shape), dtype=self.radii.dtype)
            dataset[:] = self.radii[:]
            f.create_group("temperature_increments")

    def WriteResults(self, filename, process_info):
        step = process_info[KratosMultiphysics.STEP]
        time = step = process_info[KratosMultiphysics.TIME]
        # Open the HDF5 file.
        with h5py.File(filename, "a") as f:
            assert self.radii.shape == self.temperature_increments.shape
            # Create a dataset to store the radii and temperatures data.
            dataset = f["/temperature_increments"].create_dataset(
                str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype
            )
            # Write the radii and temperatures data to the dataset. TODO: the radii is not being written, I believe
            dataset[:] = self.temperature_increments
            # Add a time label to the dataset.
            dataset.attrs["time"] = time

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

    # TODO: move out of the solver class into a sort of utility class or file?
    def TableToFunction(self, table):
        """
        Takes a Kratos PiecewiseLinearTable and outputs a Python function that accesses it.

        Parameters
        ----------
        table: Kratos PiecewiseLinearTable

        Returns
        -------
        TableAsFunction: Python function
        """

        def TableAsFunction(x):
            return table.GetValue(x)

        return TableAsFunction

    def ExtendFunction(self, f, xmax, xmin=0):
        """
        Extends a function f defined on the interval (xmin, xmax) by setting its value to zero outside
        this interval. The returned function f_extended is a piecewise function that returns f(x) for
        xmin <= x <= xmax and 0 for x < xmin or x > xmax.

        Parameters
        ----------
        f : callable
            The function to extend, defined on the open interval (xmin, xmax). Must take a single
            float input and return a float output.
        xmin : float
            The lower bound of the interval where f is defined.
        xmax : float
            The upper bound of the interval where f is defined, must be greater than xmin.

        Returns
        -------
        callable
            The extended function f_extended, defined for all real numbers. Returns f(x) for
            xmin <= x <= xmax and 0 otherwise.

        Raises
        ------
        ValueError
            If f is not callable, xmin >= xmax, or f(xmin) or f(xmax) are non-finite.
        TypeError
            If xmin or xmax are not numbers.
        """
        # Validate inputs
        if not callable(f):
            raise ValueError("f must be callable")
        if not isinstance(xmin, (int, float)) or not isinstance(xmax, (int, float)):
            raise TypeError("xmin and xmax must be numbers")
        if not np.isfinite(xmin) or not np.isfinite(xmax):
            raise ValueError("xmin and xmax must be finite")
        if xmin >= xmax:
            raise ValueError("xmin must be less than xmax")

        # Test function evaluation at boundaries
        try:
            f_xmin = f(xmin)
            f_xmax = f(xmax)
        except Exception as e:
            raise ValueError(f"Function evaluation failed at boundaries: {str(e)}")

        if not np.isfinite(f_xmin) or not np.isfinite(f_xmax):
            raise ValueError("f(xmin) and f(xmax) must be finite")

        # Define the extended function
        def f_extended(x):
            """
            Piecewise function that extends f(x) to zero outside [xmin, xmax].

            Args
                x : float
                    Input value to evaluate.

            Returns
                float
                    f(x) if xmin <= x <= xmax, 0 otherwise.
            """
            if xmin <= x <= xmax:
                return f(x)
            return 0.0

        return f_extended

    # TODO: move out of the solver class into a sort of utility class or file?
    def NormalizeAxisymmetricFunction(self, f, rmin=0, rmax=np.inf, numerical_zero=1e-15):
        """
        Normalizes a function f=f(r) representing the radial part of an axisymmetric function F(r, phi) = f(r)
        with 0 <= r < +inf,  0 < phi < 2pi.
        It generates a function g = g(r) such that the integral of G(r, phi) = g(r) over the real plane
        equals one.

        Note: f itself is not normalized, F is:
        Integral(F, real plane) = Integral(F(r,phi) r dr dphi, 0 < r < +inf, 0 < phi < 2 pi)
            = 2pi * Integral(f(r) r dr, 0 < r < +inf)

        Parameters
        ----------
        f: callable
            The input function to normalize, taking a single numeric input.

        Returns
        -------
        g: callable
            A function g such that Integral(g, 0 < r < +inf, 0 < phi < 2 pi) = 1.

        Raises
        ------
        ValueError: If the integral of r*f(r) over (0, +inf) is zero or infinite.
        """

        # Compute the integral of r * f(r) over (0, +inf)
        radial_integral, _ = quad(lambda r: r * f(r), rmin, rmax)

        # Check if the integral is finite and non-zero
        if not np.isfinite(radial_integral):
            raise ValueError("Integral of r*f(r) over (0, +inf) is infinite or NaN")
        if abs(radial_integral) < numerical_zero:  # Small threshold to avoid division by zero
            raise ValueError("Integral of r*f(r) over (0, +inf) is zero or nearly zero")

        normalization_factor = 2 * np.pi * radial_integral

        # Define the normalized function
        def g(r):
            return f(r) / normalization_factor

        return g

    def LoadFluenceFromTable(self, fluence_table_filename, fluence_table_max_r, normalize=False):
        """
        Load a laser pulse fluence profile from a CSV file containing (x, y) pairs, convert it to a function,
        and extend it to return zero beyond the specified maximum radius. Optionally normalize the profile
        to represent the total energy of the laser pulse.

        The CSV file should contain at least two columns, where x typically represents radial position (r)
        and y represents the radial part of fluence.

        The function is extended to return 0 for r > fluence_table_max_r.

        If normalize=True, the function is scaled to integrate to the total pulse energy (self.Q).

        Parameters
        ----------
        fluence_table_filename : str
            Path to the CSV file containing (x, y) pairs for the fluence profile.
        fluence_table_max_r : float
            Maximum radial position beyond which the fluence is set to zero.
        normalize : bool, optional
            If True, normalize the fluence profile to integrate to 1 and scale by self.Q (pulse energy).
            If False, use the raw fluence values from the table. Default is False.

        Returns
        -------
        FluenceFunction: callable
            A function that takes a position dictionary (with key 'r') and parameters, returning the
            fluence value at that radial position.

        Raises
        ------
        AttributeError
            If self.Q is not defined when normalize=True.
        RuntimeError
            If the table reading or function conversion fails in Kratos utilities.

        Notes
        -----
        - Assumes the CSV file uses comma as delimiter and has numeric data.
        - The returned function expects position as a dictionary with key 'r' (radial position).
        """

        # Configure table settings
        fluence_table_settings = KratosMultiphysics.Parameters("""{
            "name"             : "csv_table",
            "filename"         : "",
            "delimiter"        : ",",
            "skiprows"         : 0,
            "first_column_id"  : 0,
            "second_column_id" : 1,
            "table_id"         : -1,
            "na_replace"       : 0.0
        }""")
        fluence_table_settings["filename"].SetString(fluence_table_filename)

        # Read table
        try:
            FluenceTable = ReadCsvTableUtility(fluence_table_settings).Read()
        except Exception as e:
            raise RuntimeError(f"Failed to read fluence table: {str(e)}")

        # Convert to function
        FluenceTableAsFunction = self.TableToFunction(FluenceTable)
        if not callable(FluenceTableAsFunction):
            raise RuntimeError("TableToFunction did not return a callable function")

        # Extend function
        FluenceTableExtended = self.ExtendFunction(FluenceTableAsFunction, xmax=fluence_table_max_r)
        if not callable(FluenceTableExtended):
            raise RuntimeError("ExtendFunction did not return a callable function")

        if normalize:
            # Check for pulse energy
            if not hasattr(self, "Q") or not isinstance(self.Q, (int, float)) or self.Q <= 0:
                raise AttributeError("self.Q must be a positive scalar when normalize=True")

            # Normalize
            FluenceFunctionNormalized = self.NormalizeAxisymmetricFunction(
                FluenceTableExtended, rmin=0, rmax=fluence_table_max_r
            )
            if not callable(FluenceFunctionNormalized):
                raise RuntimeError("NormalizeAxisymmetricFunction did not return a callable function")

            def FluenceFunction(position, parameters):
                return self.Q * FluenceFunctionNormalized(position["r"])
        else:

            def FluenceFunction(position, parameters):
                return FluenceTableExtended(position["r"])

        return FluenceFunction
