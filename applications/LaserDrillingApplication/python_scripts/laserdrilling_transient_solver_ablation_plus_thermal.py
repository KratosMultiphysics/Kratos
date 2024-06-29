import numpy as np

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

    def ResidualHeatStage(self):
        pass

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
        self.ImposeTemperatureIncreaseDueTo1DConduction()
        super().RemoveElementsByAblation()
