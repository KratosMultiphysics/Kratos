import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator

from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model = cs_tools.CreateInstance(self.parameters["settings"]["model"])
        self.omega = self.settings["omega"].GetDouble()

    def Initialize(self):
        super().Initialize()

        self.model.size_in = self.model.size_out = self.x.GetNumpyArray().shape[0]
        self.model.out = self.x.deepcopy()
        self.model.Initialize()
        self.components += [self.model]

    def SolveSolutionStep(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        xt = self.solver_wrappers[1].SolveSolutionStep(y)
        r = xt - self.x
        self.model.Add(r, xt)
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            if not self.model.IsReady():
                dx = self.omega * r
            else:
                dx = self.model.Predict(r) + r
            self.x += dx
            y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(y)
            r = xt - self.x
            self.model.Add(r, xt)
            self.FinalizeIteration(r)
        # self.model.pn = y.GetNumpyArray() / 1000

