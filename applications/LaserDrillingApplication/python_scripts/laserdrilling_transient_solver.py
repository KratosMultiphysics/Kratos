import os
import numpy as np
from scipy.interpolate import interp1d
import h5py
import time as timer

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
        self.jump_between_pulses_counter = 0
        self.starting_time = timer.time()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        '''original_delta_time = self.ComputeDeltaTime()
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        time_jump_between_pulses = self.time_jump_between_pulses
        adaptative_delta_time = original_delta_time * (0.5 + np.abs(np.sin(np.pi * current_time / (2.0 * time_jump_between_pulses))))
        print("\nAdaptative_delta_time:", adaptative_delta_time, "\n")
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, adaptative_delta_time)'''
        #self.settings["time_stepping"]["time_step"].SetDouble(2e-6)
        #print("Delta time after:", self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME])

        self.jump_between_pulses_counter += delta_time
        if self.jump_between_pulses_counter >= self.time_jump_between_pulses:
            self.jump_between_pulses_counter = 0
            self.pulse_number += 1
            print("\nPulse_number:", self.pulse_number)
            if self.ablation_energy_fraction:
                self.RemoveElementsByAblation()
            if self.residual_heat_fraction:
                self.ImposeTemperatureDistributionDueToLaser1D()

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
            elem.Set(KratosMultiphysics.ACTIVE, True)

        for node in self.main_model_part.Nodes:
            node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 0.0)

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
            self.residual_heat_fraction = 0.05
        else:
            self.residual_heat_fraction = project_parameters["problem_data"]["residual_heat_fraction"].GetDouble()
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

        self.R_far = mask_aperture_diameter * 0.5
        self.cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        self.kappa = self.conductivity / (self.rho * self.cp)

        self.plot_decomposed_volume_graph = False
        self.element_id_to_study = 1578
        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []

        self.V = 4.72e-7 # mm3. Approximate ablated volume for 1 pulses (experimental). For 5 pulses it should be around 2.36e-6

        self.ablation_energy_fraction = 1 - self.residual_heat_fraction
        self.sigma = 0.5 * self.R_far
        self.pulse_number = 0
        self.K = 1 / (2 * self.sigma**2)
        self.C = self.ablation_energy_fraction * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
        A = np.pi * self.R_far**2

        if self.ablation_energy_fraction:
            # Find F_th_fraction multiplying F_th so Radius_th = R_far
            # This gives a F_th_fraction of 0.313, a little too flat (maximum not captured)
            # F_th_fraction = self.C * np.exp(-self.K * self.R_far**2) / (self.ablation_energy_fraction * self.Q / A)

            F_th_fraction = 0.333

            self.F_th = F_th_fraction * self.ablation_energy_fraction * self.Q / A
            F_th_in_cm_units = 100 * self.F_th
            print("\nF_th in J/cm2:", F_th_in_cm_units)
            self.radius_th = np.sqrt(np.log(self.C / self.F_th) / self.K)
            self.l_s = self.PenetrationDepth() # In mm
            print("\nl_s in mm:", self.l_s)
            l_s_in_nm = 1e6 * self.l_s  # In nm
            print("\nl_s in nm:", l_s_in_nm)

            # self.l_s = 0.002

        if os.path.exists("temperature_alpha.txt"):
            os.remove("temperature_alpha.txt")
        self.temperature_alpha_file = open("temperature_alpha.txt", "a")
        if os.path.exists("time_alpha.txt"):
            os.remove("time_alpha.txt")
        self.time_alpha_file = open("time_alpha.txt", "a")
        if os.path.exists("decomposed_volume_evolution.txt"):
            os.remove("decomposed_volume_evolution.txt")
        self.decomposed_volume_file = open("decomposed_volume_evolution.txt", "a")

        # TODO: Initial condition, ambient temperature. Do this using GUI!
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.T0)
            #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, self.T0)

        initial_system_energy = self.MonitorEnergy()
        #print("\nInitial system energy:", initial_system_energy)

        self.pulse_number += 1
        print("\nPulse_number:", self.pulse_number)
        if self.ablation_energy_fraction:
            self.RemoveElementsByAblation()

        computed_energy_after_laser = self.MonitorEnergy()
        #print("\nComputed energy after laser:", computed_energy_after_laser)

        residual_heat = self.residual_heat_fraction * self.Q
        print("\nResidual_heat:", residual_heat)

        if self.residual_heat_fraction:
            self.ImposeTemperatureDistributionDueToLaser1D()
        
        system_energy_after_residual_heat = self.MonitorEnergy()
        #print("System energy after residual heat:", system_energy_after_residual_heat)

        difference_after_and_before_residual_heat = system_energy_after_residual_heat - computed_energy_after_laser
        print("Difference after and before residual heat:", difference_after_and_before_residual_heat, "\n")

        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2
        self.near_field_nodes = [node for node in self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def ImposeTemperatureDistributionDueToLaser1D(self):
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        self.minimum_characteristic_Z =  1e6
        self.maximum_characteristic_Z = -1e6
        for node in self.main_model_part.Nodes:
            radius = node.Y
            if not self.ablation_energy_fraction:
                Z_interp = 0.0
            else:
                F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
                Z_interp = F(radius)
            z = node.X - Z_interp
            if radius <= self.R_far:
                delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
                old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                new_temp = old_temp + delta_temp
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)
        print("\nResidual heat fraction:", self.residual_heat_fraction)
        print("\nMinimum characteristic depth:", self.minimum_characteristic_Z)
        print("\nMaximum characteristic depth:", self.maximum_characteristic_Z)
        problem_characteristic_time_minimum_depth = self.minimum_characteristic_Z**2 / self.kappa
        problem_characteristic_time_maximum_depth = self.maximum_characteristic_Z**2 / self.kappa
        minimum_time_step_for_minimum_depth = 0.1 * problem_characteristic_time_minimum_depth
        minimum_time_step_for_maximum_depth = 0.1 * problem_characteristic_time_maximum_depth
        print("\nThermal problem characteristic time for minimum depth:", problem_characteristic_time_minimum_depth)
        print("\nThermal problem characteristic time for maximum depth:", problem_characteristic_time_maximum_depth)
        print("\nNecessary time step for minimum depth:", minimum_time_step_for_minimum_depth)
        print("\nNecessary time step for maximum depth:", minimum_time_step_for_maximum_depth, '\n')

    def EnergyPerUnitArea1D(self, radius):
        C = (1 - self.ablation_energy_fraction) * self.Q * self.K / (np.pi * (1 - np.exp(-self.K * self.R_far**2)))
        q =  C * np.exp(-self.K * radius**2)
        return q

    def TimeToAchieveEvaporationTemperature(self, radius):
        # 4.0 (2**2) in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        C = 4.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * (self.T_e - self.T0)**2)
        t = C * self.EnergyPerUnitArea1D(radius)**2
        return t

    def TemperatureVariationInZDueToLaser1D(self, radius, z):
        q = self.EnergyPerUnitArea1D(radius)
        t_evap = self.TimeToAchieveEvaporationTemperature(radius)
        # 2.0 in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        C = 2.0 * q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_evap))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_evap))
        characteristic_Z = np.sqrt(4.0 * self.kappa * t_evap)
        if characteristic_Z <= self.minimum_characteristic_Z:
            self.minimum_characteristic_Z = characteristic_Z
        if characteristic_Z >= self.maximum_characteristic_Z:
            self.maximum_characteristic_Z = characteristic_Z
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
        elapsed_time = timer.time() - self.starting_time
        print("\nElapsed_time:", elapsed_time, '\n')

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
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        for elem in self.main_model_part.Elements:
            X_centroid = elem.GetGeometry().Center().X
            Y_centroid = elem.GetGeometry().Center().Y
            if self.pulse_number == 1:
                X_interp = 0
            else:
                F = interp1d(Y, X, bounds_error=False)
                X_interp = F(Y_centroid)
            DeltaX = X_centroid - X_interp
            d_ev = self.EvaporationDepth(Y_centroid)
            if DeltaX <= d_ev and Y_centroid <= self.radius_th:
                elem.Set(KratosMultiphysics.ACTIVE, False)
                for node in elem.GetNodes():
                    node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)
        for elem in self.main_model_part.Elements:
            number_of_decomposed_nodes = 0
            for node in elem.GetNodes():
                if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                    number_of_decomposed_nodes +=1
            if number_of_decomposed_nodes == 3:
                elem.Set(KratosMultiphysics.ACTIVE, False)
        self.AddDecomposedNodesToSurfaceList()
        print('\nR_far:', self.R_far)
        print('\nRadius_th:', self.radius_th)
        print("\nDecomposed volume:", self.MonitorDecomposedVolume())

    def sortSecond(self, val):
        return val[1]

    def AddDecomposedNodesToSurfaceList(self):
        self.list_of_decomposed_nodes_ids = []
        self.list_of_decomposed_nodes_coords = []
        self.list_of_decomposed_nodes_coords_X = []
        self.list_of_decomposed_nodes_coords_Y = []
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                        self.list_of_decomposed_nodes_ids.append(node.Id)
        self.list_of_decomposed_nodes_ids = list(set(self.list_of_decomposed_nodes_ids))

        for node in self.main_model_part.Nodes:
            if node.Id in self.list_of_decomposed_nodes_ids:
                X = node.X
                Y = node.Y
                coords = [X, Y]
                self.list_of_decomposed_nodes_coords.append(coords)
        self.list_of_decomposed_nodes_coords.sort(key=self.sortSecond)
        self.list_of_decomposed_nodes_coords_X = [coord[0] for coord in self.list_of_decomposed_nodes_coords]
        self.list_of_decomposed_nodes_coords_Y = [coord[1] for coord in self.list_of_decomposed_nodes_coords]
        if os.path.exists("list_of_decomposed_nodes_coords.txt"):
            os.remove("list_of_decomposed_nodes_coords.txt")
        self.decomposed_nodes_coords_file = open("list_of_decomposed_nodes_coords.txt", "a")
        for coord in self.list_of_decomposed_nodes_coords:
            self.decomposed_nodes_coords_file.write(str(coord[0]) + " " + str(coord[1]) + "\n")
        self.decomposed_nodes_coords_file.close()

    def PenetrationDepth(self):
        F_th = self.F_th
        V = self.V
        R_th = self.radius_th
        l_s = V / (np.pi * (0.5 * R_th**2 * (np.log(self.C / F_th)) - 0.25 * self.K * R_th**4))
        return l_s

    def EvaporationDepth(self, r):
        if r >= self.radius_th:
            return 0
        else:
            q = self.C * np.exp(-self.K * r**2)
            d_ev = 0.5 * self.l_s * (np.log(q) - np.log(self.F_th))
            return d_ev

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
