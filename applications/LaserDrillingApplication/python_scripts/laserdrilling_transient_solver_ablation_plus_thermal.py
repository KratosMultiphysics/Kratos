import numpy as np

from KratosMultiphysics.LaserDrillingApplication import laserdrilling_transient_solver

def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolverAblationPlusThermal(model, custom_settings)

class LaserDrillingTransientSolverAblationPlusThermal(laserdrilling_transient_solver.LaserDrillingTransientSolver):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

    def SolveSolutionStep(self):
        #print("\nLaserDrillingTransientSolverAblationPlusThermal strategy!!!\n")
        super(laserdrilling_transient_solver.LaserDrillingTransientSolver, self).SolveSolutionStep()

    def ResidualHeatStage(self):
        pass

    def TemperatureVariationInZDueToLaser1D(self, radius, z):
        q = self.EnergyPerUnitArea1D(radius)
        t_penetration = self.InitialThermalConductionTime(radius)
        # 2.0 in numerator due to equation (4c) in Weber, 2014. 'Heat accumulation during pulsed laser materials processing'
        C = 2.0 * q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_penetration))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_penetration))
        return delta_temp

    def RemoveElementsByAblation(self):
        self.ImposeTemperatureIncreaseDueTo1DConduction()
        super().RemoveElementsByAblation()
