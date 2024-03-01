import os
import numpy as np
import h5py

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver

def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolver(model, custom_settings)

class LaserDrillingTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    """The transient class for laser convection-diffusion solvers.

    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_transient_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the settings in base class
        super().__init__(model, custom_settings)
        self.jump_between_pulses_counter = 0

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.jump_between_pulses_counter += delta_time
        if self.jump_between_pulses_counter >= self.time_jump_between_pulses:
            self.jump_between_pulses_counter = 0
            initial_computed_energy = self.MonitorEnergy()
            print()
            print("initial computed energy:", initial_computed_energy)
            self.ImposeTemperatureDistributionDueToLaser1D()
            computed_energy_after_laser = self.MonitorEnergy()
            print()
            print("computed energy after laser:", computed_energy_after_laser)
            print()
            print("laser energy:", self.Q)
            diff_of_energies = computed_energy_after_laser - initial_computed_energy
            print()
            print("diff between final and initial energies:", diff_of_energies)
            print()
            computed_energy_after_adjusting = self.MonitorEnergy()
            target_total_energy = initial_computed_energy + self.Q
            print()
            print("target total energy:", target_total_energy)
            print()
            print("computed energy after adjusting:", computed_energy_after_adjusting)
            relative_error = 100 * (computed_energy_after_adjusting - target_total_energy) / target_total_energy
            print()
            print("relative error:", relative_error, "%")
            print()

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

    def Initialize(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).Initialize()

        # Set element counter variable to zero
        for elem in self.main_model_part.Elements:
            elem.SetValue(LaserDrillingApplication.THERMAL_COUNTER, 0)
            elem.SetValue(KratosMultiphysics.TEMPERATURE, 0.0)
            elem.SetValue(LaserDrillingApplication.THERMAL_DECOMPOSITION, 0.01)
            elem.SetValue(LaserDrillingApplication.DECOMPOSED_ELEMENTAL_VOLUME, 0.0)

        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        material_settings = materials["properties"][0]["Material"]

        with open("ProjectParameters.json", 'r') as project_parameters_file:
            project_parameters = KratosMultiphysics.Parameters(project_parameters_file.read())

        if not project_parameters["problem_data"].Has("energy"):
            self.Q = 25e-6
        else:
            self.Q = project_parameters["problem_data"]["energy"].GetDouble()
        if not project_parameters["problem_data"].Has("residual_heat_fraction"):
            residual_heat_fraction = 0.05
        else:
            residual_heat_fraction = project_parameters["problem_data"]["residual_heat_fraction"].GetDouble()
        if not project_parameters["problem_data"].Has("mask_aperture_diameter"):
            mask_aperture_diameter = 0.025
        else:
            mask_aperture_diameter = project_parameters["problem_data"]["mask_aperture_diameter"].GetDouble()
        if not project_parameters["problem_data"].Has("vaporisation_temperature"):
            self.T_e = 693.0
        else:
            self.T_e = project_parameters["problem_data"]["vaporisation_temperature"].GetDouble()
        if not project_parameters["problem_data"].Has("time_jump_between_pulses"):
            self.time_jump_between_pulses = 1e6
        else:
            self.time_jump_between_pulses = project_parameters["problem_data"]["time_jump_between_pulses"].GetDouble()

        self.Q *= residual_heat_fraction
        self.R_far = mask_aperture_diameter * 0.5
        self.cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        self.kappa = self.conductivity / (self.rho * self.cp)
        self.plot_decomposed_volume_graph = False

        self.element_id_to_study = 1804

        '''print(project_parameters["problem_data"]["end_time"].GetDouble())
        project_parameters["problem_data"]["end_time"].SetDouble(2e-5)
        print(project_parameters["problem_data"]["end_time"].GetDouble()
        stop'''

        # TODO: Initial condition, ambient temperature. Do this using GUI!
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
            #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)

        initial_computed_energy = self.MonitorEnergy()
        print()
        print("initial computed energy:", initial_computed_energy)

        self.ImposeTemperatureDistributionDueToLaser1D()

        computed_energy_after_laser = self.MonitorEnergy()
        print()
        print("computed energy after laser:", computed_energy_after_laser)

        print()
        print("laser energy:", self.Q)
        diff_of_energies = computed_energy_after_laser - initial_computed_energy
        print()
        print("diff between final and initial energies:", diff_of_energies)

        '''temperature_factor = self.Q / diff_of_energies
        print()
        print("temperature_factor:", temperature_factor)
        self.AdjustTemperatureField(temperature_factor)'''

        computed_energy_after_adjusting = self.MonitorEnergy()
        target_total_energy = initial_computed_energy + self.Q
        print()
        print("target total energy:", target_total_energy)
        print()
        print("computed energy after adjusting:", computed_energy_after_adjusting)
        print()
        print("actual laser energy added:", computed_energy_after_adjusting - initial_computed_energy)
        print()
        relative_error = 100 * (computed_energy_after_adjusting - target_total_energy) / target_total_energy
        print()
        print("relative error:", relative_error, "%")
        print()

        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2
        self.near_field_nodes = [node for node in self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

        import os
        if os.path.exists("decomposed_volume_evolution.txt"):
            os.remove("decomposed_volume_evolution.txt")
        self.decomposed_volume_file = open("decomposed_volume_evolution.txt", "a")
        if os.path.exists("temperature_alpha.txt"):
            os.remove("temperature_alpha.txt")
        self.temperature_alpha_file = open("temperature_alpha.txt", "a")
        if os.path.exists("time_alpha.txt"):
            os.remove("time_alpha.txt")
        self.time_alpha_file = open("time_alpha.txt", "a")

    def ImposeTemperatureDistributionDueToLaser1D(self):
        for node in self.main_model_part.Nodes:
            radius = node.Y
            z = node.X
            if radius <= self.R_far:
                delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
                old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temp = old_temp + delta_temp
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)

    def AdjustTemperatureField(self, temperature_factor):
        for node in self.main_model_part.Nodes:
            radius = node.Y
            if radius <= self.R_far:
                temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                temp *= temperature_factor
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

    def EnergyPerUnitArea1D(self, radius):

        sigma = 0.5 * self.R_far
        C = self.Q / (2 * np.pi * sigma**2 * (1 - np.exp(-self.R_far**2 / (2 * sigma**2))))       
        q =  C * np.exp(-0.5 * (radius / self.R_far)**2)
        return q

    def TimeToAchieveEvaporationTemperature(self, radius):

        C = 1.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * (self.T_e - self.T0)**2)
        t = C * self.EnergyPerUnitArea1D(radius)**2
        return t

    def TemperatureVariationInZDueToLaser1D(self, radius, z):

        q = self.EnergyPerUnitArea1D(radius)
        t_evap = self.TimeToAchieveEvaporationTemperature(radius)
        C = q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_evap))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_evap))
        return delta_temp

    def SolveSolutionStep(self):
        super(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver, self).SolveSolutionStep()
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

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

    def Finalize(self):
        super().Finalize()
        self.decomposed_volume_file.close()
        self.temperature_alpha_file.close()
        self.time_alpha_file.close()

        if self.plot_decomposed_volume_graph:
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
            #Y0 = elem.GetNode(0).Y
            #Y1 = elem.GetNode(1).Y
            #Y2 = elem.GetNode(2).Y
            #if Y0 <= self.R_far or Y1 <= self.R_far or Y2 <= self.R_far:
            out = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.THERMAL_ENERGY, self.main_model_part.ProcessInfo)
            # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
            energy += out[0]

        return energy

    def ImposeTemperatureDueToLaser(self):

        t_ini = 5e-3
        self.C_L = 2.0 * self.Q / (8.0 * self.cp * (np.pi * self.kappa)**1.5 * self.rho)

        def bell_curve(t, radius_squared):
            z = -radius_squared / (4.0 * self.kappa * t)
            bell_curve_value = self.C_L / t**1.5 * np.exp(z)
            return bell_curve_value

        for node in self.main_model_part.Nodes:
            t_ini = 5e-3
            r_2 = node.X**2 + node.Y**2 + node.Z**2
            temp = self.T0 + bell_curve(t_ini, r_2)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

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
            assert self.radii.shape  == self.temperature_increments.shape

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
