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

    def SetParameters(self):
        super().SetParameters()

        ## 2024 Woodfield - Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue
        ## Laser data
        self.omega_0 = self.R_far # mm
        self.F_p = 2.0 * self.Q / (np.pi * self.omega_0**2) # J/mm2

        ## Material calibration using experiments
        if not self.material_settings["Variables"].Has("OPTICAL_PENETRATION_DEPTH"):
            self.delta_pen = 5e-4 # mm
        else:
            self.delta_pen = self.material_settings['Variables']['OPTICAL_PENETRATION_DEPTH'].GetDouble()
        if not self.material_settings["Variables"].Has("IONIZATION_ENERGY_PER_VOLUME_THRESHOLD"):
            self.q_ast = 10.0 # J/mm3
        else:
            self.q_ast = self.material_settings['Variables']['IONIZATION_ENERGY_PER_VOLUME_THRESHOLD'].GetDouble()

        self.minimum_of_two_energies_per_volume = self.q_ast

        ##
        #self.delta_pen = self.l_s
        self.r_ast_max = self.omega_0 * np.sqrt(0.5 * np.log(self.F_p / (self.delta_pen * self.q_ast)))

        # Compute ionization energy per volume of C11H12O3
        E_m_H = 1312e3 #  J/mol (1st level ionization energy)
        E_m_C = 4621e3 #  J/mol (3rd level ionization energy)
        E_m_O = 3388e3 #  J/mol (2nd level ionization energy)
        W_m   = 192e-3 # Kg/mol (molecular weight)
        self.ionization_energy_per_volume_threshold = self.rho * (E_m_H + E_m_C + E_m_O) / W_m # J/mm3

    def ImposeTemperatureIncreaseDueTo1DConduction(self):
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        for node in self.main_model_part.Nodes:
            radius = node.Y
            F = interp1d(Y, X, bounds_error=False, fill_value=0.0)
            distance_to_surface = F(radius)
            z = node.X - distance_to_surface
            delta_temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
            old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            new_temp = old_temp + delta_temp
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)

            delta_pen = self.delta_pen
            F_p = self.F_p
            omega_0 = self.omega_0
            q_energy_per_volume = (1.0 / delta_pen) * F_p * np.exp(- 2.0 * (radius / omega_0)**2) * np.exp(- z / delta_pen)
            node.SetValue(LaserDrillingApplication.ENERGY_PER_VOLUME, q_energy_per_volume)

            # Compute enthalpy energy per volume
            delta_temp = self.T_e - old_temp
            enthalpy_energy_per_volume = self.rho * (self.H_ev + self.cp * delta_temp)
            node.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume)

        for elem in self.main_model_part.Elements:
            q_energy_per_volume = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo)
            elem.SetValue(LaserDrillingApplication.ENERGY_PER_VOLUME, q_energy_per_volume[0])
            enthalpy_energy_per_volume = elem.CalculateOnIntegrationPoints(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo)
            elem.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume[0])

    def TemperatureVariationInZDueToLaser1D(self, radius, z):        
        delta_pen = self.delta_pen
        F_p = self.F_p
        omega_0 = self.omega_0
        q_energy_per_volume = (1.0 / delta_pen) * F_p * np.exp(- 2.0 * (radius / omega_0)**2) * np.exp(- z / delta_pen)
        delta_temp = q_energy_per_volume / (self.rho * self.cp)            
        return delta_temp

    def EvaporationDepth(self, r):
        if r >= self.r_ast_max:
            return 0.0
        else:
            delta_pen = self.delta_pen
            F_p = self.F_p
            q_ast = self.q_ast
            omega_0 = self.omega_0
            z_ast = delta_pen * (np.log(F_p / (delta_pen * q_ast)) - 2.0 * (r / omega_0)**2)
            return z_ast

    def RemoveElementsByAblation(self):

        '''initial_system_energy = self.MonitorEnergy()
        print("\n\nEnergy before laser:", initial_system_energy)'''

        self.ImposeTemperatureIncreaseDueTo1DConduction()

        '''expected_energy_after_laser = initial_system_energy + self.Q
        initial_system_energy = self.MonitorEnergy()
        print("Expected energy after laser:", expected_energy_after_laser)
        print("Actual energy after laser:", initial_system_energy)'''

        self.RemoveElementsUsingEnergyPerVolumeThreshold()

        '''initial_system_energy = self.MonitorEnergy()
        print("Energy after ablation:", initial_system_energy)
        decomp_vol = self.MonitorDecomposedVolume()
        print("Actual volume loss due to laser:",  decomp_vol)

        delta_pen = self.delta_pen
        F_p = self.F_p
        q_ast = self.minimum_of_two_energies_per_volume #self.q_ast
        omega_0 = self.omega_0
        vol_n_pulses = self.pulse_number * 0.25 * delta_pen * np.pi * omega_0**2 * (np.log(F_p / (delta_pen * q_ast)))**2
        print("Expected volume loss due to laser:", vol_n_pulses, "\n\n")'''

    def ResidualHeatStage(self):
        pass

    def RemoveElementsUsingEnergyPerVolumeThreshold(self):
        if self.ablation_energy_fraction:
            for elem in self.main_model_part.Elements:
                q_energy_per_volume = elem.GetValue(LaserDrillingApplication.ENERGY_PER_VOLUME)
                enthalpy_energy_per_volume = elem.GetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME)
                ionization_energy_per_volume_threshold = self.ionization_energy_per_volume_threshold
                energy_threshold = min(enthalpy_energy_per_volume, ionization_energy_per_volume_threshold)
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
