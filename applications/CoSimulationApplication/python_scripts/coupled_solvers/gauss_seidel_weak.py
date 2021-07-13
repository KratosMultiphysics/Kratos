# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(settings, models, solver_name):
    return GaussSeidelWeakCoupledSolver(settings, models, solver_name)

class GaussSeidelWeakCoupledSolver(CoSimulationCoupledSolver):
    def SolveSolutionStep(self):
        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        for solver_name, solver in self.solver_wrappers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()
            # solving only the coupling operations of the current solver
            for coupling_op in self.coupling_operations_dict.values():
                if coupling_op.settings["solver"].GetString() == solver_name:
                    coupling_op.FinalizeSolveSolutionStep()
            self._SynchronizeOutputData(solver_name)

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        return True
