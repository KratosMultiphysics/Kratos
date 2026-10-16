# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(settings, models, solver_name):
    return GaussSeidelWeakCoupledSolver(settings, models, solver_name)

class GaussSeidelWeakCoupledSolver(CoSimulationCoupledSolver):
    def SolveSolutionStep(self):
        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        for data_transfer_operator in self.data_transfer_operators_dict.values():
            data_transfer_operator.InitializeNonLinearIteration()

        for solver_name, solver in self.solver_wrappers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver_name)

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        for data_transfer_operator in self.data_transfer_operators_dict.values():
            data_transfer_operator.FinalizeNonLinearIteration()

        return True
