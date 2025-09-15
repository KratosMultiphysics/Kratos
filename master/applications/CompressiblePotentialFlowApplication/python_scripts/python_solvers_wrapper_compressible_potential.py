import KratosMultiphysics

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_type == "potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_solver
        return potential_flow_solver.CreateSolver(model, solver_settings)
    elif solver_type == "ale_potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import ale_potential_flow_solver
        return ale_potential_flow_solver.CreateSolver(model, solver_settings, parallelism)
    elif solver_type == "adjoint_potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_adjoint_solver
        return potential_flow_adjoint_solver.CreateSolver(model, solver_settings)
    else:
        err_msg = "The requested solver type is not in the python solvers wrapper. Solver type is : " + solver_type
        raise ValueError(err_msg)


def CreateSolver(model, custom_settings):
    if not isinstance(model, KratosMultiphysics.Model):
        raise TypeError("input is expected to be provided as a Kratos Model object")

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise TypeError("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
