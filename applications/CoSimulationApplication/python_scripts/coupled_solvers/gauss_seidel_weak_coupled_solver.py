from __future__ import print_function, absolute_import, division

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(model, cosim_solver_settings, solver_name):
    return GaussSeidelWeakCouplingSolver(model, cosim_solver_settings, solver_name)

class GaussSeidelWeakCouplingSolver(CoSimulationCoupledSolver):
    def SolveSolutionStep(self):
        for coupling_op in self.coupling_operations_list:
            coupling_op.InitializeCouplingIteration()

        for solver_name, solver in self.participating_solvers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver_name)

        for coupling_op in self.coupling_operations_list:
            coupling_op.FinalizeCouplingIteration()

        return True
