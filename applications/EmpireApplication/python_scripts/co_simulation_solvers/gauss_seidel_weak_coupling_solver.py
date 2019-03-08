from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver

# Other imports
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import couplingsolverprint, green

def CreateSolver(cosim_solver_settings, level):
    return GaussSeidelWeakCouplingSolver(cosim_solver_settings, level)

class GaussSeidelWeakCouplingSolver(CoSimulationBaseCouplingSolver):
    def SolveSolutionStep(self):
        for solver_name in self.solver_names:
            solver = self.solvers[solver_name]
            self._SynchronizeInputData(solver, solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver, solver_name)

    def _Name(self):
        return self.__class__.__name__
