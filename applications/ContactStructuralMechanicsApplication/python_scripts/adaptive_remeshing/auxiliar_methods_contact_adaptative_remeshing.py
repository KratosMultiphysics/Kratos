# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import auxiliar methods
from KratosMultiphysics.StructuralMechanicsApplication.adaptive_remeshing.auxiliar_methods_structural_adaptative_remeshing import AuxiliarMethodsStructuralAdaptiveRemeshing

class AuxiliarMethodsContactAdaptiveRemeshing(AuxiliarMethodsStructuralAdaptiveRemeshing):
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

    def SPRAdaptativeRemeshingRunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered with SPR convergence criteria

            Keyword arguments:
            self It signifies an instance of a class.
        """

        # Remeshing adaptively
        solver = self.analysis._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        convergence_criteria = solver.get_convergence_criterion()
        builder_and_solver = solver.get_builder_and_solver()
        mechanical_solution_strategy = solver.get_mechanical_solution_strategy()

        while self.analysis.KeepAdvancingSolutionLoop():
            self.analysis.time = solver.AdvanceInTime(self.analysis.time)
            non_linear_iteration = 1
            while non_linear_iteration <= self.analysis.non_linear_iterations:
                if non_linear_iteration == 1 or self.analysis._CheckIfModelIsModified():
                    self.analysis.InitializeSolutionStep()
                    solver.Predict()
                    self.analysis._ResetModelIsModified()
                    self.analysis.is_printing_rank = False
                computing_model_part.ProcessInfo.SetValue(KM.NL_ITERATION_NUMBER, non_linear_iteration)
                is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                solver.SolveSolutionStep()
                is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
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
                    # Compute metric
                    metric_process = solver.get_metric_process()
                    metric_process.Execute()

                    # Clean up contact pairs
                    CSMA.ContactUtilities.CleanContactModelParts(computing_model_part)

                    # We remove the contact submodelparts
                    computing_model_part.RemoveSubModelPart("Contact")

                    # Ensure properties defined
                    KM.AuxiliarModelPartUtilities(computing_model_part.GetRootModelPart()).RecursiveEnsureModelPartOwnsProperties(True)

                    # We create the contact submodelparts
                    computing_model_part.CreateSubModelPart("Contact")

                    # Remeshing
                    remeshing_process = solver.get_remeshing_process()
                    remeshing_process.Execute()

                    self.analysis._SetModelIsModified()
                    non_linear_iteration += 1
            self.analysis.OutputSolutionStep()
