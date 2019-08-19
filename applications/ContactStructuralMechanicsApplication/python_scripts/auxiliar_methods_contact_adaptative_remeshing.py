from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.MeshingApplication as MA

# Import auxiliar methods
from KratosMultiphysics.auxiliar_methods_adaptative_remeshing import AuxiliarMethodsAdaptiveRemeshing

class AuxiliarMethodsContactAdaptiveRemeshing(AuxiliarMethodsAdaptiveRemeshing):
    """
    This class is an auxiliar script when using adaptative remeshing put in a class for contact

    It can be imported and used as "black-box"
    """
    def __init__(self, analysis):
        """This function is the constructor of the class

            Keyword arguments:
            self It signifies an instance of a class.
            analysis The AnalysisStage to be computed
        """
        super(AuxiliarMethodsContactAdaptiveRemeshing, self).__init__(analysis)

    def ExecuteBeforeFinalizeSolutionStep(self):
        """This function is executed before the FinalizeSolutionStep

            Keyword arguments:
            self It signifies an instance of a class.
        """
        solver = self.analysis._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        map_parameters = KM.Parameters("""
        {
            "echo_level"                       : 0,
            "origin_variable_historical"       : false,
            "destination_variable_historical"  : false,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "origin_variable"                  : "AUGMENTED_NORMAL_CONTACT_PRESSURE",
            "destination_variable"             : "AUGMENTED_NORMAL_CONTACT_PRESSURE",
            "integration_order"                : 2
        }
        """)

        interface_model_part = computing_model_part.GetSubModelPart("Contact")
        main_model_part = computing_model_part.GetRootModelPart()
        main_model_part.RemoveSubModelPart("ComputingContact")

        if interface_model_part.HasSubModelPart("SlaveSubModelPart"):
            slave_interface_model_part = interface_model_part.GetSubModelPart("SlaveSubModelPart")
        else:
            slave_interface_model_part = interface_model_part.CreateSubModelPart("SlaveSubModelPart")
            KM.FastTransferBetweenModelPartsProcess(slave_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES, KM.SLAVE).Execute()
            KM.FastTransferBetweenModelPartsProcess(slave_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS, KM.SLAVE).Execute()
        if interface_model_part.HasSubModelPart("MasterSubModelPart"):
            master_interface_model_part = interface_model_part.GetSubModelPart("MasterSubModelPart")
        else:
            master_interface_model_part = interface_model_part.CreateSubModelPart("MasterSubModelPart")
            KM.FastTransferBetweenModelPartsProcess(master_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES, KM.MASTER).Execute()
            KM.FastTransferBetweenModelPartsProcess(master_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS, KM.MASTER).Execute()

        mortar_mapping = KM.SimpleMortarMapperProcess(slave_interface_model_part, master_interface_model_part, map_parameters)
        mortar_mapping.Execute()

    def SPRAdaptativeRemeshingRunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered with SPR convergence criteria

            Keyword arguments:
            self It signifies an instance of a class.
        """

        # Remeshing adaptively
        solver = self.analysis._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        root_model_part = computing_model_part.GetRootModelPart()
        convergence_criteria = solver.get_convergence_criterion()
        builder_and_solver = solver.get_builder_and_solver()
        mechanical_solution_strategy = solver.get_mechanical_solution_strategy()

        while self.analysis.KeepAdvancingSolutionLoop():
            self.analysis.time = solver.AdvanceInTime(self.analysis.time)
            non_linear_iteration = 1
            while non_linear_iteration <= self.analysis.non_linear_iterations:
                if root_model_part.Is(KM.MODIFIED):
                    self.analysis.ClearDatabase()
                    self.analysis.ReInitializeSolver()
                if non_linear_iteration == 1 or root_model_part.Is(KM.MODIFIED):
                    self.analysis.InitializeSolutionStep()
                    solver.Predict()
                    root_model_part.Set(KM.MODIFIED, False)
                    self.analysis.is_printing_rank = False
                computing_model_part.ProcessInfo.SetValue(KM.NL_ITERATION_NUMBER, non_linear_iteration)
                is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                solver.SolveSolutionStep()
                is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                self.ExecuteBeforeFinalizeSolutionStep()
                self.analysis.FinalizeSolutionStep()
                if is_converged:
                    self.analysis.is_printing_rank = True
                    KM.Logger.PrintInfo(self.analysis._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                    break
                elif non_linear_iteration == self.analysis.non_linear_iterations:
                    self.analysis.is_printing_rank = True
                    KM.Logger.PrintInfo(self.analysis._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                    break
                else:
                    # Before remesh we set the flag INTERFACE to the conditions (we need edges to preserve submodelparts)
                    KM.VariableUtils().SetFlag(KM.INTERFACE, True, computing_model_part.GetSubModelPart("Contact").Conditions)

                    # We remove the contact model part to avoid problems (it will  be recomputed later)
                    contact_model_part = computing_model_part.GetSubModelPart("Contact")
                    for model_part in contact_model_part.SubModelParts:
                        contact_model_part.RemoveSubModelPart(model_part.Name)
                    computing_model_part.RemoveSubModelPart("ComputingContact")

                    # Ensure properties
                    MA.MeshingUtilities.EnsureModelPartOwnsProperties(root_model_part)
                    MA.MeshingUtilities.EnsureModelPartOwnsProperties(computing_model_part)

                    metric_process = solver.get_metric_process()
                    remeshing_process = solver.get_remeshing_process()
                    metric_process.Execute()
                    remeshing_process.Execute()

                    # We remove the contact model part to avoid problems (it will  be recomputed later)
                    computing_model_part.RemoveSubModelPart("Contact")

                    root_model_part.Set(KM.MODIFIED, True)
                    non_linear_iteration += 1
            self.analysis.OutputSolutionStep()
