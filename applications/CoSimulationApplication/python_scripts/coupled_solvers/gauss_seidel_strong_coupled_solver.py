from __future__ import print_function, absolute_import, division

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCouplingSolver

# Other imports
from KratosMultiphysics.CoSimulationApplication.convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import couplingsolverprint, red, green, cyan, bold
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(model, cosim_solver_settings):
    return GaussSeidelStrongCouplingSolver(model, cosim_solver_settings)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, model, cosim_solver_settings):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(GaussSeidelStrongCouplingSolver, self).__init__(model, cosim_solver_settings)

        self.convergence_accelerators_list = cs_tools.CreateConvergenceAccelerators(
            self.cosim_solver_settings["convergence_accelerators"],
            self.solvers)
        # self.convergence_accelerator.SetEchoLevel(self.echo_level) # TODO set echo-lvl?

        self.convergence_criteria = CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria"],
            self.solvers)
        self.convergence_criteria.SetEchoLevel(self.echo_level)

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]

    def Initialize(self):
        super(GaussSeidelStrongCouplingSolver, self).Initialize()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Initialize()
        self.convergence_criteria.Initialize()

    def Finalize(self):
        super(GaussSeidelStrongCouplingSolver, self).Finalize()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Finalize()
        self.convergence_criteria.Finalize()

    def InitializeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).InitializeSolutionStep()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.InitializeSolutionStep()
        self.convergence_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).FinalizeSolutionStep()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.FinalizeSolutionStep()
        self.convergence_criteria.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.InitializeNonLinearIteration()
            self.convergence_criteria.InitializeNonLinearIteration()

            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.FinalizeNonLinearIteration()
            self.convergence_criteria.FinalizeNonLinearIteration()

            if self.convergence_criteria.IsConverged():
                if self.echo_level > 0:
                    couplingsolverprint(self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                break
            else:
                for conv_acc in self.convergence_accelerators_list:
                    conv_acc.ComputeUpdate()

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))

    def PrintInfo(self):
        super(GaussSeidelStrongCouplingSolver, self).PrintInfo()

        couplingsolverprint(self._Name(), "Uses the following objects:")
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.PrintInfo()
        self.convergence_criteria.PrintInfo()

    def Check(self):
        super(GaussSeidelStrongCouplingSolver, self).Check()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Check()
        self.convergence_criteria.Check()

    def _Name(self):
        return self.__class__.__name__