from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import time
import json

cs_data_structure = cs_tools.cs_data_structure
import matplotlib.pyplot as plt

def Create(parameters):
    return CoupledSolverRelaxation(parameters)


class CoupledSolverRelaxation(CoupledSolverGaussSeidel):
    def SolveSolutionStep(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        xt = self.solver_wrappers[1].SolveSolutionStep(y)
        r = xt - self.x
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            self.x += self.omega * r
            y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(y)
            r = xt - self.x
            self.FinalizeIteration(r)
