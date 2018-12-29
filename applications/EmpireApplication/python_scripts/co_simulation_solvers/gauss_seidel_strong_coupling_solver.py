from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver

# Other imports
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import couplingsolverprint, red, green, cyan, bold

def CreateSolver(cosim_solver_settings, level):
    return GaussSeidelStrongCouplingSolver(cosim_solver_settings, level)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, cosim_solver_settings, level):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(GaussSeidelStrongCouplingSolver, self).__init__(cosim_solver_settings, level)

        self.convergence_accelerator = CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"],
            self.solvers, self.lvl)
        self.convergence_accelerator.SetEchoLevel(self.echo_level)

        self.convergence_criteria = CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"],
            self.solvers, self.lvl)
        self.convergence_criteria.SetEchoLevel(self.echo_level)

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]

    def Initialize(self):
        super(GaussSeidelStrongCouplingSolver, self).Initialize()
        self.convergence_accelerator.Initialize()
        self.convergence_criteria.Initialize()

    def Finalize(self):
        super(GaussSeidelStrongCouplingSolver, self).Finalize()
        self.convergence_accelerator.Finalize()
        self.convergence_criteria.Finalize()

    def InitializeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).InitializeSolutionStep()
        self.convergence_accelerator.InitializeSolutionStep()
        self.convergence_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).FinalizeSolutionStep()
        self.convergence_accelerator.FinalizeSolutionStep()
        self.convergence_criteria.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            self.convergence_accelerator.InitializeNonLinearIteration()
            self.convergence_criteria.InitializeNonLinearIteration()

            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)

            self.convergence_accelerator.FinalizeNonLinearIteration()
            self.convergence_criteria.FinalizeNonLinearIteration()

            if self.convergence_criteria.IsConverged():
                if self.echo_level > 0:
                    couplingsolverprint(self.lvl, self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                break
            else:
                self.convergence_accelerator.ComputeUpdate()

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))

    def PrintInfo(self):
        super(GaussSeidelStrongCouplingSolver, self).PrintInfo()

        couplingsolverprint(self.lvl, self._Name(), "Uses the following objects:")
        self.convergence_accelerator.PrintInfo()
        self.convergence_criteria.PrintInfo()

    def Check(self):
        super(GaussSeidelStrongCouplingSolver, self).Check()
        self.convergence_accelerator.Check()
        self.convergence_criteria.Check()

    def _Name(self):
        return self.__class__.__name__