from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

def AdaptativeRemeshingRunSolutionLoop(analysis):
    """This function executes the solution loop of the AnalysisStage for cases where remeshing may be considered

        Keyword arguments:
        analysis The AnalysisStage to be reinitialized
    """

    # If we remesh using a process
    computing_model_part = analysis._GetSolver().GetComputingModelPart()
    root_model_part = computing_model_part.GetRootModelPart()

    while analysis.time < analysis.end_time:
        analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
        # We reinitialize if remeshed previously
        if root_model_part.Is(KratosMultiphysics.MODIFIED):
            ReInitializeSolver(analysis)
        analysis.InitializeSolutionStep()
        # We reinitialize if remeshed on the InitializeSolutionStep
        if root_model_part.Is(KratosMultiphysics.MODIFIED):
            ReInitializeSolver(analysis)
            analysis.InitializeSolutionStep()
        analysis._GetSolver().Predict()
        analysis._GetSolver().SolveSolutionStep()
        analysis.FinalizeSolutionStep()
        analysis.OutputSolutionStep()

def ReInitializeSolver(analysis):
    """ This reinitializes after remesh

        Keyword arguments:
        analysis The AnalysisStage to be reinitialized
    """
    analysis._GetSolver().Clear()
    # WE INITIALIZE THE SOLVER
    analysis._GetSolver().Initialize()
    # WE RECOMPUTE THE PROCESSES AGAIN
    ## Processes initialization
    for process in analysis._GetListOfProcesses():
        process.ExecuteInitialize()
    ## Processes before the loop
    for process in analysis._GetListOfProcesses():
        process.ExecuteBeforeSolutionLoop()
    ## Processes of initialize the solution step
    for process in analysis._GetListOfProcesses():
        process.ExecuteInitializeSolutionStep()
