
import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Dynamic" or solver_type == "dynamic"):
            time_integration_method = solver_settings["time_integration_method"].GetString()
            if (time_integration_method == "implicit"):
                solver_module_name = "mpm_implicit_dynamic_solver"
            elif (time_integration_method == "explicit"):
                solver_module_name = "mpm_explicit_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\", \"explicit\""
                raise Exception(err_msg)
        elif (solver_type == "Quasi-static" or solver_type == "quasi_static"):
            solver_module_name = "mpm_quasi_static_solver"
        elif (solver_type == "Static" or solver_type == "static"):
            solver_module_name = "mpm_static_solver"
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\", \"quasi_static\""
            raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\""
        raise Exception(err_msg)

    # Remove settings that are not needed any more
    solver_settings.RemoveValue("solver_type")
    solver_settings.RemoveValue("time_integration_method") # does not throw even if the value is not existing

    module_full = 'KratosMultiphysics.MPMApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver


def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)