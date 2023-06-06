import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model: KratosMultiphysics.Model, solver_settings: KratosMultiphysics.Parameters, parallelism: str):

    filter_type = solver_settings["filter_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):

        if (filter_type == "general_scalar" ):
            solver_module_name = "helmholtz_scalar_solver"

        elif (filter_type == "bulk_surface_shape" or filter_type == "general_vector"):
            solver_module_name = "helmholtz_vector_solver"
        else:
            err_msg =  'The requested solver type "' + filter_type + '" is not in the python solvers wrapper\n'
            err_msg += 'Available options are: "general_scalar", "general_vector", "bulk_surface_shape"'
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):

        if (filter_type == filter_type == "bulk_surface_shape" or filter_type == "general_vector"):
            solver_module_name = "trilinos_general_vector_filter_solver"

        elif (filter_type == "general_scalar" ):
            solver_module_name = "trilinos_general_scalar_filter_solver"

        else:
            err_msg =  'The requested solver type "' + filter_type + '" is not in the python solvers wrapper\n'
            err_msg += 'Available options are: "trilinos_general_vector_filter_solver", "trilinos_general_scalar_filter_solver", "trilinos_shape_filter_solver"'
            raise Exception(err_msg)

    else:
        err_msg =  'The requested parallel type "' + parallelism + '" is not available!\n'
        err_msg += 'Available options are: "OpenMP", "MPI"'
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.OptimizationApplication.filtering.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model: KratosMultiphysics.Model, custom_settings: KratosMultiphysics.Parameters):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
