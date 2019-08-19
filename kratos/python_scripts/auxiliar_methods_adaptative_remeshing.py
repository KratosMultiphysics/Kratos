from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class AuxiliarMethodsAdaptiveRemeshing(object):
    """
    This class is an auxiliar script when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, analysis):
        """This function is the constructor of the class

            Keyword arguments:
            self It signifies an instance of a class.
            analysis The AnalysisStage to be computed
        """
        # Saves the analysis
        self.analysis = analysis

    def AdaptativeRemeshingDetectBoundary(self):
        """This function detects the boundary to preserve pure node BC

            Keyword arguments:
            self It signifies an instance of a class.
        """
        solver = self.analysis._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        if not self.analysis.process_remesh:
            convergence_criteria = solver.get_convergence_criterion()
            convergence_criteria.Initialize(computing_model_part)

        # Ensuring to have conditions on the BC before remesh
        is_surface = False
        for elem in computing_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            break

        if not is_surface:
            list_model_parts = []
            # We need to detect the conditions in the boundary conditions
            if self.analysis.project_parameters.Has("constraints_process_list"):
                constraints_process_list = self.analysis.project_parameters["constraints_process_list"]
                for i in range(0,constraints_process_list.size()):
                    item = constraints_process_list[i]
                    list_model_parts.append(item["Parameters"]["model_part_name"].GetString())
            skin_detection_parameters = KratosMultiphysics.Parameters("""
            {
                "list_model_parts_to_assign_conditions" : []
            }
            """)
            root_model_part_name = computing_model_part.GetRootModelPart().Name
            for name in list_model_parts:
                name = name.replace(root_model_part_name + ".", "")
                skin_detection_parameters["list_model_parts_to_assign_conditions"].Append(name)

            if computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                detect_skin = KratosMultiphysics.SkinDetectionProcess2D(computing_model_part, skin_detection_parameters)
            else:
                detect_skin = KratosMultiphysics.SkinDetectionProcess3D(computing_model_part, skin_detection_parameters)
            detect_skin.Execute()
        solver.SetEchoLevel(self.analysis.echo_level)

    def AdaptativeRemeshingRunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered

            Keyword arguments:
            self It signifies an instance of a class.
        """

        # If we remesh using a process
        computing_model_part = self.analysis._GetSolver().GetComputingModelPart()
        root_model_part = computing_model_part.GetRootModelPart()

        while self.analysis.KeepAdvancingSolutionLoop():
            self.analysis.time = self.analysis._GetSolver().AdvanceInTime(self.analysis.time)
            # We reinitialize if remeshed previously
            if root_model_part.Is(KratosMultiphysics.MODIFIED):
                self.analysis.ReInitializeSolver()
            self.analysis.InitializeSolutionStep()
            # We reinitialize if remeshed on the InitializeSolutionStep
            if root_model_part.Is(KratosMultiphysics.MODIFIED):
                self.analysis.ReInitializeSolver()
                self.analysis.InitializeSolutionStep()
            self.analysis._GetSolver().Predict()
            self.analysis._GetSolver().SolveSolutionStep()
            self.analysis.FinalizeSolutionStep()
            self.analysis.OutputSolutionStep()

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
                if root_model_part.Is(KratosMultiphysics.MODIFIED):
                    self.analysis.ReInitializeSolver()
                if non_linear_iteration == 1 or root_model_part.Is(KratosMultiphysics.MODIFIED):
                    self.analysis.InitializeSolutionStep()
                    solver.Predict()
                    root_model_part.Set(KratosMultiphysics.MODIFIED, False)
                    self.analysis.is_printing_rank = False
                computing_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, non_linear_iteration)
                is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                solver.SolveSolutionStep()
                is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                self.ExecuteBeforeFinalizeSolutionStep()
                self.analysis.FinalizeSolutionStep()
                if is_converged:
                    self.analysis.is_printing_rank = True
                    KratosMultiphysics.Logger.PrintInfo(self.analysis._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                    break
                elif non_linear_iteration == self.analysis.non_linear_iterations:
                    self.analysis.is_printing_rank = True
                    KratosMultiphysics.Logger.PrintInfo(self.analysis._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                    break
                else:
                    # Remesh
                    metric_process = solver.get_metric_process()
                    remeshing_process = solver.get_remeshing_process()
                    metric_process.Execute()
                    remeshing_process.Execute()

                    root_model_part.Set(KratosMultiphysics.MODIFIED, True)
                    non_linear_iteration += 1
            self.analysis.OutputSolutionStep()

    def ExecuteBeforeFinalizeSolutionStep(self):
        """This function is executed before the FinalizeSolutionStep

            Keyword arguments:
            self It signifies an instance of a class.
        """
        pass
