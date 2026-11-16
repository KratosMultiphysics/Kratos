import KratosMultiphysics


def ValidateMechanicalSettings(mechanical_settings):
    """Reject mechanical configurations incompatible with the coupling.

    Only a quasi-static Newton-Raphson mechanical solve is supported. The check
    is shared by the solvers and by the analysis stage so that an unsupported
    configuration is reported before the strategy is constructed.
    """
    solution_type = mechanical_settings["solution_type"].GetString()
    if solution_type != "Quasi-Static":
        raise RuntimeError(
            "DamApplication: hydro-mechanical uplift coupling is only "
            "supported for 'solution_type' == 'Quasi-Static' "
            "(got '{}').".format(solution_type))

    strategy_type = mechanical_settings["strategy_type"].GetString()
    if strategy_type != "Newton-Raphson":
        raise RuntimeError(
            "DamApplication: hydro-mechanical uplift coupling requires "
            "'mechanical_solver_settings.strategy_type' == 'Newton-Raphson' "
            "(got '{}').".format(strategy_type))


def SolveMechanicalStepWithUpliftCoupling(
        mechanical_strategy, update_uplift, coupling_settings):
    """Solve one time step with the uplift pressure coupled to the joint state.

    The mechanical problem is solved to convergence, the uplift pressure is
    recomputed from the resulting joint aperture and the two are iterated until
    the pressure, the aperture and the crack front stop changing. The
    constitutive state of the step is committed only once, after convergence.

    ``update_uplift(iteration)`` must update the pressure loads and return
    ``(is_converged, metrics)``.
    """
    max_iterations = coupling_settings["max_iterations"].GetInt()
    echo_level = coupling_settings["echo_level"].GetInt()

    if max_iterations < 1:
        raise RuntimeError(
            "DamApplication: 'hydro_mechanical_coupling.max_iterations' must be "
            "at least 1.")

    mechanical_strategy.InitializeSolutionStep()
    mechanical_strategy.Predict()

    metrics = None
    for iteration in range(1, max_iterations + 1):
        mechanical_strategy.SolveSolutionStep()

        is_converged, metrics = update_uplift(iteration)

        if echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo(
                "DamUpliftCoupling",
                "iteration {:d}: {}".format(iteration, metrics))

        if is_converged:
            mechanical_strategy.FinalizeSolutionStep()
            if echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo(
                    "DamUpliftCoupling",
                    "converged after {} iteration(s): {}".format(
                        iteration, metrics))
            return iteration, metrics

    raise RuntimeError(
        "DamApplication: the hydro-mechanical uplift coupling did not converge "
        "in {} iterations. Last state: {}. The time step is NOT advanced with "
        "an inconsistent pressure field.".format(max_iterations, metrics))
