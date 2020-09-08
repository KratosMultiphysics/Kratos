# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.factories.helpers as factories_helper
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(settings, models, solver_name):
    return GaussSeidelStrongCoupledSolver(settings, models, solver_name)

class GaussSeidelStrongCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        self.convergence_accelerators_list = factories_helper.CreateConvergenceAccelerators(
            self.settings["convergence_accelerators"],
            self.solver_wrappers,
            self.echo_level)

        self.convergence_criteria_list = factories_helper.CreateConvergenceCriteria(
            self.settings["convergence_criteria"],
            self.solver_wrappers,
            self.echo_level)

        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt()

    def Initialize(self):
        super().Initialize()

        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Initialize()

        for conv_crit in self.convergence_criteria_list:
            conv_crit.Initialize()

    def Finalize(self):
        super().Finalize()

        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Finalize()

        for conv_crit in self.convergence_criteria_list:
            conv_crit.Finalize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for conv_acc in self.convergence_accelerators_list:
            conv_acc.InitializeSolutionStep()

        for conv_crit in self.convergence_criteria_list:
            conv_crit.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        for conv_acc in self.convergence_accelerators_list:
            conv_acc.FinalizeSolutionStep()

        for conv_crit in self.convergence_criteria_list:
            conv_crit.FinalizeSolutionStep()


    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), colors.cyan("Coupling iteration:"), colors.bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            for coupling_op in self.coupling_operations_dict.values():
                coupling_op.InitializeCouplingIteration()

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.InitializeNonLinearIteration()

            for conv_crit in self.convergence_criteria_list:
                conv_crit.InitializeNonLinearIteration()

            for solver_name, solver in self.solver_wrappers.items():
                self._SynchronizeInputData(solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)

            for coupling_op in self.coupling_operations_dict.values():
                coupling_op.FinalizeCouplingIteration()

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.FinalizeNonLinearIteration()

            for conv_crit in self.convergence_criteria_list:
                conv_crit.FinalizeNonLinearIteration()

            is_converged = all([conv_crit.IsConverged() for conv_crit in self.convergence_criteria_list])

            if is_converged:
                if self.echo_level > 0:
                    cs_tools.cs_print_info(self._ClassName(), colors.green("### CONVERGENCE WAS ACHIEVED ###"))
                self.__CommunicateStateOfConvergence(True)
                return True

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), colors.red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))
                self.__CommunicateStateOfConvergence(True) # True because max number of iterations is achieved. Otherwise external solver is stuck in time
                return False

            # if it reaches here it means that the coupling has not converged and this was not the last coupling iteration
            self.__CommunicateStateOfConvergence(False)

            # do relaxation only if this iteration is not the last iteration of this timestep
            for conv_acc in self.convergence_accelerators_list:
                conv_acc.ComputeAndApplyUpdate()


    def Check(self):
        super().Check()

        if len(self.convergence_criteria_list) == 0:
            raise Exception("At least one convergence criteria has to be specified")

        # TODO check if an accelerator was specified for a field that is manipulated in the input!

        for conv_crit in self.convergence_criteria_list:
            conv_crit.Check()

        for conv_crit in self.convergence_accelerators_list:
            conv_crit.Check()

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "convergence_accelerators" : [],
            "convergence_criteria"     : [],
            "num_coupling_iterations"  : 10
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultSettings())

        return this_defaults

    def __CommunicateStateOfConvergence(self, is_converged):
        # Communicate the state of convergence with external solvers through IO
        convergence_signal_config = {
            "type" : "convergence_signal",
            "is_converged" : is_converged
        }

        for solver in self.solver_wrappers.values():
            solver.ExportData(convergence_signal_config)

