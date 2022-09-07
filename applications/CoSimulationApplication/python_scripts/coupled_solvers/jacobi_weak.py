# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(settings, models, solver_name, base_class = CoSimulationCoupledSolver): 
    return GetJacobiWeakCoupledSolver(base_class)(settings, models, solver_name)

def GetJacobiWeakCoupledSolver(base):
    class TemplatedJacobiWeakCoupledSolver(base):
        def SolveSolutionStep(self):
            for coupling_op in self.coupling_operations_dict.values():
                coupling_op.InitializeCouplingIteration()

            for solver_name, solver in self.solver_wrappers.items():
                self._SynchronizeInputData(solver_name)

            for solver_name, solver in self.solver_wrappers.items():
                solver.SolveSolutionStep()

            for solver_name, solver in self.solver_wrappers.items():
                self._SynchronizeOutputData(solver_name)

            for coupling_op in self.coupling_operations_dict.values():
                coupling_op.FinalizeCouplingIteration()

            return True

    return TemplatedJacobiWeakCoupledSolver