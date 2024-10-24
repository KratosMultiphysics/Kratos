# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel_strong import GaussSeidelStrongCoupledSolver

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(settings, models, solver_name):
    return BlockStrongCoupledSolver(settings, models, solver_name)

class BlockStrongCoupledSolver(GaussSeidelStrongCoupledSolver):
    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            self.process_info[KratosCoSim.COUPLING_ITERATION_NUMBER] += 1

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

                # Apply relaxation for each solver output
                for conv_acc in self.convergence_accelerators_list:
                    conv_acc.ComputeAndApplyUpdate()

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
                self.__CommunicateIfTimeStepNeedsToBeRepeated(False)
                return True

            if k+1 >= self.num_coupling_iterations:
                if self.echo_level > 0:
                    cs_tools.cs_print_info(self._ClassName(), colors.red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))
                self.__CommunicateIfTimeStepNeedsToBeRepeated(False)
                return False

            # if it reaches here it means that the coupling has not converged and this was not the last coupling iteration
            self.__CommunicateIfTimeStepNeedsToBeRepeated(True)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "convergence_accelerators" : [],
            "convergence_criteria"     : [],
            "num_coupling_iterations"  : 10
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

        return this_defaults

    def __CommunicateIfTimeStepNeedsToBeRepeated(self, repeat_time_step):
        # Communicate if the time step needs to be repeated with external solvers through IO
        export_config = {
            "type" : "repeat_time_step",
            "repeat_time_step" : repeat_time_step
        }

        for solver in self.solver_wrappers.values():
            solver.ExportData(export_config)
