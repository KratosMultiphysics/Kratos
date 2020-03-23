from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

cs_data_structure = cs_tools.cs_data_structure

def Create(parameters):
    return CoupledSolverRelaxation(parameters)


class CoupledSolverRelaxation(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.omega = self.settings["omega"].GetDouble()

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
