import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_type == "potential_flow":
        import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver
        return potential_flow_solver.CreateSolver(model, solver_settings)
    elif solver_type == "ale_potential_flow":
        import KratosMultiphysics.CompressiblePotentialFlowApplication.ale_potential_flow_solver as ale_potential_flow_solver
        return ale_potential_flow_solver.CreateSolver(model, solver_settings, parallel_type)
    elif solver_type == "adjoint_potential_flow":
        import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_adjoint_solver as adjoint_solver
        return adjoint_solver.CreateSolver(model, solver_settings)
    else:
        raise Exception("Solver type '"+str(solver_type)+"' not added. Please specify an available solver")

    return solver
