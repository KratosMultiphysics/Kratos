def Factory(*_args, **_kwargs):
    raise RuntimeError(
        "\"optimization_problem_restart_input_process\" can no longer be configured under "
        "\"processes\"; restart is now configured via the top-level \"restart_settings\" section "
        "of the optimization parameters (see OptimizationAnalysis.GetDefaultParameters()). "
        "Restore itself is no longer a process at all -- it happens directly in "
        "OptimizationAnalysis.__init__()."
    )
