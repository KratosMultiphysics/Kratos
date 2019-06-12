from __future__ import print_function, absolute_import, division

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCouplingSolver

# Other imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import couplingsolverprint, red, green, cyan, bold
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(model, cosim_solver_settings, solver_name):
    return GaussSeidelStrongCouplingSolver(model, cosim_solver_settings, solver_name)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, model, cosim_solver_settings, solver_name):
        if not cosim_solver_settings['coupling_sequence'].size() == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(GaussSeidelStrongCouplingSolver, self).__init__(model, cosim_solver_settings, solver_name)

        self.convergence_accelerators_list = cs_tools.CreateConvergenceAccelerators(
            self.settings["convergence_accelerators"],
            self.participating_solvers,
            self.echo_level)

        self.convergence_criteria_list = cs_tools.CreateConvergenceCriteria(
            self.settings["convergence_criteria"],
            self.participating_solvers,
            self.echo_level)

        self.coupling_operations_list.extend(self.convergence_accelerators_list)
        self.coupling_operations_list.extend(self.convergence_criteria_list)

        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt()


    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self._Name(), cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            for coupling_op in self.coupling_operations_list:
                coupling_op.InitializeCouplingIteration()

            for solver_name, solver in self.participating_solvers.items():
                self._SynchronizeInputData(solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)

            for coupling_op in self.coupling_operations_list:
                coupling_op.FinalizeCouplingIteration()

            is_converged = True
            for conv_crit in self.convergence_criteria_list:
                is_converged = is_converged and conv_crit.IsConverged()
            if is_converged:
                if self.echo_level > 0:
                    couplingsolverprint(self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                return True
            else:
                # TODO I think this should not be done in the last iterations if the solution does not converge in this timestep
                for conv_acc in self.convergence_accelerators_list:
                    conv_acc.Execute()

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))
                return False

    def Check(self):
        # TODO check if at least one conv-crit was specified?
        # TODO check if an accelerator was specified for a field that is manipulated in the input!
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = cs_tools.cs_data_structure.Parameters("""{
            "convergence_accelerators" : [],
            "convergence_criteria"     : [],
            "num_coupling_iterations"  : 10
        }""")
        this_defaults.AddMissingParameters(super(GaussSeidelStrongCouplingSolver, cls)._GetDefaultSettings())

        return this_defaults
