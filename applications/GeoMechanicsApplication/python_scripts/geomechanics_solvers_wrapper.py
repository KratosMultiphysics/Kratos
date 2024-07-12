import KratosMultiphysics
from importlib import import_module


def CreateSolverByParameters(model, custom_settings, parallelism):
    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type_raw = custom_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        solver_type = solver_type_raw.lower()
        if solver_type in ("u_pw", "geomechanics_u_pw_solver", "twophase"):
            solver_module_name = "geomechanics_U_Pw_solver"

        elif solver_type in ("pw", "geomechanics_pw_solver"):
            solver_module_name = "geomechanics_Pw_solver"

        elif solver_type in ("t", "geomechanics_t_solver"):
            solver_module_name = "geomechanics_T_solver"
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"geomechanics_U_Pw_solver\", \"geomechanics_Pw_solver\", \"geomechanics_T_solver\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\""
        raise Exception(err_msg)

    module_full_name = 'KratosMultiphysics.GeoMechanicsApplication.' + solver_module_name
    solver = import_module(module_full_name).CreateSolver(model, custom_settings)
    return solver


def CreateSolver(model, custom_settings):
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    return CreateSolverByParameters(model, custom_settings["solver_settings"], parallelism)