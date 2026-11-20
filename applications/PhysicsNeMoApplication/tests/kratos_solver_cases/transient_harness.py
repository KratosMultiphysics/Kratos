"""Shared time-loop harness for the in-memory transient solver cases.

AnalysisStage.RunSolutionLoop's own sequence, opened up so a state can be
collected after each converged step (and so the loop body can be
instrumented) - the real time loop the app's tests previously faked with
ProcessInfo[STEP] assignments.
"""

import numpy


def RunTransientAnalysis(analysis, collect=None):
    """Runs a Kratos analysis step by step, collecting a state per step.

    The real time loop the app's tests previously faked with
    ProcessInfo[STEP] assignments - AnalysisStage.RunSolutionLoop's own
    sequence, opened up so a state can be collected after each converged
    step (and so the loop body can be instrumented).

    Args:
        analysis: An AnalysisStage instance (not yet Run/Initialized).
        collect: callable(model_part) -> array-like state, or None to
            collect nothing.

    Returns:
        (T, ...) numpy array of the collected states (empty array when
        collect is None).
    """
    analysis.Initialize()
    model_part = analysis._GetSolver().GetComputingModelPart()

    states = []
    while analysis.KeepAdvancingSolutionLoop():
        analysis.time = analysis._AdvanceTime()
        analysis.InitializeSolutionStep()
        analysis.SolveSolutionStep()
        analysis.FinalizeSolutionStep()
        analysis.OutputSolutionStep()
        if collect is not None:
            states.append(numpy.asarray(collect(model_part), dtype=float))
    analysis.Finalize()
    return numpy.stack(states) if states else numpy.empty((0,))
