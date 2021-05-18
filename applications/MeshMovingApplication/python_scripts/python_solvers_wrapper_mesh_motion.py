import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()
    if solver_type.startswith("mesh_solver_"):
        solver_type = solver_type[12:] # remove preceding "mesh_solver_"
        depr_msg  = 'Using the old names to construct the MeshSolver\n'
        depr_msg += 'Please remove the preceding "mesh_solver_" from "solver_type"'
        KratosMultiphysics.Logger.PrintWarning("DEPRECATION", depr_msg)

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):

        if (solver_type == "laplacian"):
            solver_module_name = "mesh_solver_laplacian"

        elif (solver_type == "structural_similarity"):
            solver_module_name = "mesh_solver_structural_similarity"

        else:
            err_msg =  'The requested solver type "' + solver_type + '" is not in the python solvers wrapper\n'
            err_msg += 'Available options are: "laplacian", "structural_similarity"'
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):

        if (solver_type == "laplacian"):
            solver_module_name = "trilinos_mesh_solver_laplacian"

        elif (solver_type == "structural_similarity"):
            solver_module_name = "trilinos_mesh_solver_structural_similarity"

        else:
            err_msg =  'The requested solver type "' + solver_type + '" is not in the python solvers wrapper\n'
            err_msg += 'Available options are: "laplacian", "structural_similarity"'
            raise Exception(err_msg)

    else:
        err_msg =  'The requested parallel type "' + parallelism + '" is not available!\n'
        err_msg += 'Available options are: "OpenMP", "MPI"'
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.MeshMovingApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
