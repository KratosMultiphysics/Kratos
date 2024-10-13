import numpy as np
from scipy.interpolate import interp1d

import KratosMultiphysics
import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
from KratosMultiphysics.LaserDrillingApplication import laserdrilling_transient_solver

def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolverAblationPlusThermal(model, custom_settings)

class LaserDrillingTransientSolverAblationPlusThermal(laserdrilling_transient_solver.LaserDrillingTransientSolver):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

    def SolveSolutionStep(self):
        super(laserdrilling_transient_solver.LaserDrillingTransientSolver, self).SolveSolutionStep()

    def ComputeIonizationEnergyPerUnitVolumeThreshold(self):
        # Compute ionization energy per volume of C11_H12_O3
        E_m_H = 1312e3 #  J/mol (1st level ionization energy)
        E_m_C = 4621e3 #  J/mol (3rd level ionization energy)
        E_m_O = 3388e3 #  J/mol (2nd level ionization energy)
        W_m   = 192e-3 # Kg/mol (molecular weight)
        return self.rho * (E_m_H + E_m_C + E_m_O) / W_m # J/mm3

    def ComputeOpticalPenetrationDepth(self):
        light_lambda = 550e-6 # mm, light wavelength
        epoxy_n = self.refractive_index_n
        n = epoxy_n
        A = 4.0 * n / ((n + 1)**2 + n**2)
        self.l_s = 0.25 * light_lambda * A / np.pi
        return self.l_s

    def ImposeTemperatureIncreaseDueToLaser(self):
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        print("\nPulse number", self.pulse_number, '\n')
        for node in self.main_model_part.Nodes:
            radius = node.Y
            F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
            distance_to_surface = F(radius)
            z = node.X - distance_to_surface
            delta_temp = self.TemperatureVariationDueToLaser(radius, z)
            old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            if self.adjust_T_field_after_ablation:
                old_temp = self.reference_T_after_laser
            new_temp = old_temp + delta_temp
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)

            delta_pen = self.delta_pen
            F_p = self.F_p
            omega_0 = self.omega_0
            import math
            q_energy_per_volume = (1.0 / delta_pen) * F_p * math.exp(- 2.0 * (radius / omega_0)**2) * math.exp(- z / delta_pen)
            node.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume)

            # Compute enthalpy energy per volume
            delta_temp = self.T_e - old_temp
            enthalpy_energy_per_volume = self.rho * (self.H_ev + self.cp * delta_temp)
            node.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume)

        for elem in self.main_model_part.Elements:
            q_energy_per_volume = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo)
            elem.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume[0])
            enthalpy_energy_per_volume = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo)
            elem.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume[0])

    def TemperatureVariationDueToLaser(self, radius, z):        
        delta_pen = self.delta_pen
        F_p = self.F_p
        omega_0 = self.omega_0
        import math
        q_energy_per_volume = (1.0 / delta_pen) * F_p * math.exp(- 2.0 * (radius / omega_0)**2) * math.exp(- z / delta_pen)
        delta_temp = q_energy_per_volume / (self.rho * self.cp)            
        return delta_temp

    def RemoveElementsByAblation(self):

        '''initial_system_energy = self.MonitorEnergy()
        print("\n\nEnergy before laser:", initial_system_energy)'''

        self.ImposeTemperatureIncreaseDueToLaser()

        '''expected_energy_after_laser = initial_system_energy + self.Q
        initial_system_energy = self.MonitorEnergy()
        print("Expected energy after laser:", expected_energy_after_laser)
        print("Actual energy after laser:", initial_system_energy)'''

        self.RemoveElementsUsingEnergyPerVolumeThreshold()

        '''initial_system_energy = self.MonitorEnergy()
        print("Energy after ablation:", initial_system_energy)'''

        decomp_vol = self.MonitorDecomposedVolume()
        print("Actual volume loss due to laser:",  decomp_vol)

        delta_pen = self.delta_pen
        F_p = self.F_p
        q_ast = self.q_ast
        omega_0 = self.omega_0
        import math
        vol_n_pulses = self.pulse_number * 0.25 * delta_pen * math.pi * omega_0**2 * (math.log(F_p / (delta_pen * q_ast)))**2
        print("Expected volume loss due to laser:", vol_n_pulses, "\n")

        relative_error = 100.0 * (decomp_vol - vol_n_pulses) / vol_n_pulses
        print("Relative error in volume (%):", relative_error, "\n\n")

    def ResidualHeatStage(self):
        pass

    def RemoveElementsUsingEnergyPerVolumeThreshold(self):
        if self.ablation_energy_fraction:
            for elem in self.main_model_part.Elements:
                q_energy_per_volume = elem.GetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME)
                enthalpy_energy_per_volume = elem.GetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME)
                if self.use_enthalpy_and_ionization:
                    ionization_energy_per_volume_threshold = self.ionizarion_energy_per_volume_threshold
                    energy_threshold = min(enthalpy_energy_per_volume, ionization_energy_per_volume_threshold)
                else:
                    energy_threshold = self.q_ast
                if q_energy_per_volume >= energy_threshold:
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
    
    def Finalize(self):
        super().Finalize()
        if self.print_hole_geometry_files:
            self.hole_theoretical_profile_file = open('hole_theoretical_profile.txt', "w")
            '''for node_Y in self.list_of_decomposed_nodes_coords_Y:
                Z_coord = self.pulse_number * self.EvaporationDepth(node_Y)
                self.hole_theoretical_profile_file.write(str(node_Y) + " " + str(-Z_coord) + "\n")'''
            for i, node_Y in enumerate(self.hole_theoretical_Y_coords):
                if self.hole_theoretical_X_coords[i]:
                    self.hole_theoretical_profile_file.write(str(node_Y) + " " + str(-self.hole_theoretical_X_coords[i]) + "\n")
            self.hole_theoretical_profile_file.close()
