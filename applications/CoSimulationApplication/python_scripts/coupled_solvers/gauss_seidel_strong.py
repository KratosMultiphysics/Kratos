# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import numpy as np
import KratosMultiphysics.RomApplication as RomApp
import KratosMultiphysics.ConvectionDiffusionApplication as KCD

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

        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt()

    def Initialize(self):
        super().Initialize()

        self.convergence_accelerators_list = factories_helper.CreateConvergenceAccelerators(
            self.settings["convergence_accelerators"],
            self.solver_wrappers,
            self.data_communicator,
            self.echo_level)

        self.convergence_criteria_list = factories_helper.CreateConvergenceCriteria(
            self.settings["convergence_criteria"],
            self.solver_wrappers,
            self.data_communicator,
            self.echo_level)

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

        self.process_info[KratosCoSim.COUPLING_ITERATION_NUMBER] = 0

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        for conv_acc in self.convergence_accelerators_list:
            conv_acc.FinalizeSolutionStep()

        for conv_crit in self.convergence_criteria_list:
            conv_crit.FinalizeSolutionStep()


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
                model_part = solver._analysis_stage._GetSolver().GetComputingModelPart()
                if solver_name=="solid":
                    interface_model_part = model_part.GetSubModelPart("GENERIC_Interface_solid")
                elif solver_name=="fluid":
                    interface_model_part = model_part.GetSubModelPart("GENERIC_Interface_fluid")
                interface_nodes = set(
                    node.Id for node in interface_model_part.Nodes
                )
                self._SynchronizeInputData(solver_name)
                # print("########PRINT AFTER _SynchronizeOutputData##########")
                # if solver_name=="solid":
                #     model_part = self.solver_wrappers["solid"]._analysis_stage._GetSolver().GetComputingModelPart()
                #     interface_model_part = model_part.GetSubModelPart("GENERIC_Interface_solid")
                #     for node in interface_model_part.Nodes:
                #         face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #         print(f"Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                #         temperature = node.GetSolutionStepValue(KM.TEMPERATURE)
                #         print(f"Node ID: {node.Id}, Temperature: {np.round(temperature, 10)}")
                # for elem in model_part.Elements:
                #     if elem.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in elem.GetNodes():
                #             if node.Id in interface_nodes:
                #                 if solver_name=="fluid":
                #                     aux_flux = node.GetSolutionStepValue(KCD.AUX_FLUX)
                #                     print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Aux flux: {np.round(aux_flux, 10)}")
                #                 elif solver_name=="solid":
                #                     face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                     print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                # for cond in model_part.Conditions:
                #     if cond.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in cond.GetNodes():
                #             if node.Id in interface_nodes:
                #                 if solver_name=="fluid":
                #                     aux_flux = node.GetSolutionStepValue(KCD.AUX_FLUX)
                #                     print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Aux flux: {np.round(aux_flux, 10)}")
                #                 elif solver_name=="solid":
                #                     face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                     print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)
                # for elem in model_part.Elements:
                #     if elem.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in elem.GetNodes():
                #             if node.Id in interface_nodes:
                #                 temperature = node.GetSolutionStepValue(KM.TEMPERATURE)
                #                 print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Temperature: {np.round(temperature, 10)}")
                #                 face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                 print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                #                 reaction_flux = node.GetSolutionStepValue(KM.REACTION_FLUX)
                #                 print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Reaction flux: {np.round(reaction_flux, 10)}")
                # for cond in model_part.Conditions:
                #     if cond.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in cond.GetNodes():
                #             if node.Id in interface_nodes:
                #                 temperature = node.GetSolutionStepValue(KM.TEMPERATURE)
                #                 print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Temperature: {np.round(temperature, 10)}")
                #                 face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                 print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                #                 reaction_flux = node.GetSolutionStepValue(KM.REACTION_FLUX)
                #                 print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Reaction flux: {np.round(reaction_flux, 10)}")
                # for elem in model_part.Elements:
                #     if elem.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in elem.GetNodes():
                #             if node.Id in interface_nodes:
                #                 if solver_name=="fluid":
                #                     aux_flux = node.GetSolutionStepValue(KCD.AUX_FLUX)
                #                     print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Aux flux: {np.round(aux_flux, 10)}")
                #                 elif solver_name=="solid":
                #                     face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                     print(f"Element ID: {elem.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")
                # for cond in model_part.Conditions:
                #     if cond.GetValue(RomApp.HROM_WEIGHT)!= 1.0:
                #         for node in cond.GetNodes():
                #             if node.Id in interface_nodes:
                #                 if solver_name=="fluid":
                #                     aux_flux = node.GetSolutionStepValue(KCD.AUX_FLUX)
                #                     print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Aux flux: {np.round(aux_flux, 10)}")
                #                 elif solver_name=="solid":
                #                     face_heat_flux = node.GetSolutionStepValue(KM.FACE_HEAT_FLUX)
                #                     print(f"Condition ID: {cond.Id}, Node ID: {node.Id}, Face heat flux: {np.round(face_heat_flux, 10)}")

                print(f"################FINISHED iteration {k} ####################")


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


