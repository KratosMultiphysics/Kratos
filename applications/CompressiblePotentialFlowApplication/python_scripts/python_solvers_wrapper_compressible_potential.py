def CreateSolverByParameters(model, solver_settings, parallel_type):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_type == "potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_solver
        return potential_flow_solver.CreateSolver(model, solver_settings)
    elif solver_type == "ale_potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import ale_potential_flow_solver
        return ale_potential_flow_solver.CreateSolver(model, solver_settings, parallel_type)
    elif solver_type == "adjoint_potential_flow":
        from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_adjoint_solver
        return potential_flow_adjoint_solver.CreateSolver(model, solver_settings)
    else:
        raise ValueError("Solver type '"+str(solver_type)+"' not added. Please specify an available solver")