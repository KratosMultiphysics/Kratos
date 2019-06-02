from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

def AdaptativeRemeshingDetectBoundary(analysis):
    """This function executes the Initialize of the AnalysisStage for cases where remeshing may be considered

        Keyword arguments:
        analysis The AnalysisStage to be initialized
    """
    solver = analysis._GetSolver()
    computing_model_part = solver.GetComputingModelPart()
    if not analysis.process_remesh:
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
        if analysis.project_parameters.Has("constraints_process_list"):
            constraints_process_list = analysis.project_parameters["constraints_process_list"]
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
    solver.SetEchoLevel(analysis.echo_level)

def AdaptativeRemeshingRunSolutionLoop(analysis):
    """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered

        Keyword arguments:
        analysis The AnalysisStage to be reinitialized
    """

    # If we remesh using a process
    computing_model_part = analysis._GetSolver().GetComputingModelPart()
    root_model_part = computing_model_part.GetRootModelPart()

    while analysis.KeepAdvancingSolutionLoop():
        analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
        # We reinitialize if remeshed previously
        if root_model_part.Is(KratosMultiphysics.MODIFIED):
            analysis.ReInitializeSolver()
        analysis.InitializeSolutionStep()
        # We reinitialize if remeshed on the InitializeSolutionStep
        if root_model_part.Is(KratosMultiphysics.MODIFIED):
            analysis.ReInitializeSolver()
            analysis.InitializeSolutionStep()
        analysis._GetSolver().Predict()
        analysis._GetSolver().SolveSolutionStep()
        analysis.FinalizeSolutionStep()
        analysis.OutputSolutionStep()

def SPRAdaptativeRemeshingRunSolutionLoop(analysis, with_contact = False):
    """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered with SPR convergence criteria

        Keyword arguments:
        analysis The AnalysisStage to be reinitialized
    """

    # Remeshing adaptively
    solver = analysis._GetSolver()
    computing_model_part = solver.GetComputingModelPart()
    root_model_part = computing_model_part.GetRootModelPart()
    metric_process = solver.get_metric_process()
    remeshing_process = solver.get_remeshing_process()
    convergence_criteria = solver.get_convergence_criterion()
    builder_and_solver = solver.get_builder_and_solver()
    mechanical_solution_strategy = solver.get_mechanical_solution_strategy()

    while analysis.KeepAdvancingSolutionLoop():
        analysis.time = solver.AdvanceInTime(analysis.time)
        non_linear_iteration = 1
        while non_linear_iteration <= analysis.non_linear_iterations:
            if root_model_part.Is(KratosMultiphysics.MODIFIED):
                analysis.ReInitializeSolver()
            if non_linear_iteration == 1 or root_model_part.Is(KratosMultiphysics.MODIFIED):
                analysis.InitializeSolutionStep()
                solver.Predict()
                root_model_part.Set(KratosMultiphysics.MODIFIED, False)
                analysis.is_printing_rank = False
            computing_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, non_linear_iteration)
            is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
            solver.SolveSolutionStep()
            is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
            if with_contact:
                analysis._transfer_slave_to_master()
            analysis.FinalizeSolutionStep()
            if is_converged:
                analysis.is_printing_rank = True
                KratosMultiphysics.Logger.PrintInfo(analysis._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                break
            elif non_linear_iteration == analysis.non_linear_iterations:
                analysis.is_printing_rank = True
                KratosMultiphysics.Logger.PrintInfo(analysis._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                break
            else:
                if with_contact:
                    # Before remesh we set the flag INTERFACE to the conditions (we need edges to preserve submodelparts)
                    KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, computing_model_part.GetSubModelPart("Contact").Conditions)

                    # We remove the contact model part to avoid problems (it will  be recomputed later)
                    contact_model_part = computing_model_part.GetSubModelPart("Contact")
                    for model_part in contact_model_part.SubModelParts:
                        contact_model_part.RemoveSubModelPart(model_part.Name)
                    computing_model_part.RemoveSubModelPart("ComputingContact")

                metric_process.Execute()
                remeshing_process.Execute()

                if with_contact:
                    # We remove the contact model part to avoid problems (it will  be recomputed later)
                    computing_model_part.RemoveSubModelPart("Contact")
                root_model_part.Set(KratosMultiphysics.MODIFIED, True)
                non_linear_iteration += 1
        analysis.OutputSolutionStep()
