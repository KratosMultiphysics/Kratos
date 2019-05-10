# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver
import KratosMultiphysics
# Other imports
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import couplingsolverprint, red, green, cyan, bold


def CreateSolver(cosim_solver_settings, level):
    return SimpleSteadyCouplingSolver(cosim_solver_settings, level)

class SimpleSteadyCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, cosim_solver_settings, level):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(SimpleSteadyCouplingSolver, self).__init__(cosim_solver_settings, level)
        self.time = self.cosim_solver_settings["start_coupling_time"]

        self.convergence_criteria = CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"],
            self.solvers, self.lvl)
        self.convergence_criteria.SetEchoLevel(self.echo_level)

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]



    def SolveSolutionStep(self):

        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)
