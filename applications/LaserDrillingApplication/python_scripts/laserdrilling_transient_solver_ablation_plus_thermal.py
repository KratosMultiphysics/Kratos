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

        ##
        self.r_ast_max = self.omega_0 * np.sqrt(0.5 * np.log(self.F_p / (self.delta_pen * self.q_ast)))

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

        initial_system_energy = self.MonitorEnergy()
        print("\n\nEnergy before laser:", initial_system_energy)

        self.ImposeTemperatureIncreaseDueTo1DConduction()

        expected_energy_after_laser = initial_system_energy + self.Q
        initial_system_energy = self.MonitorEnergy()
        print("Expected energy after laser:", expected_energy_after_laser)
        print("Actual energy after laser:", initial_system_energy)

        super().RemoveElementsByAblation()

        initial_system_energy = self.MonitorEnergy()
        print("Energy after ablation:", initial_system_energy)
        decomp_vol = self.MonitorDecomposedVolume()
        print("Actual volume loss due to laser:",  decomp_vol)

        delta_pen = self.delta_pen
        F_p = self.F_p
        q_ast = self.q_ast
        omega_0 = self.omega_0
        vol_n_pulses = self.pulse_number * 0.25 * delta_pen * np.pi * omega_0**2 * (np.log(F_p / (delta_pen * q_ast)))**2
        print("Expected volume loss due to laser:", vol_n_pulses, "\n\n")

    def ResidualHeatStage(self):
        pass
