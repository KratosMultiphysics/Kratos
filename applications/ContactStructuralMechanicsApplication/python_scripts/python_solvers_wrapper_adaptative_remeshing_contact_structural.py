import KratosMultiphysics
from importlib import import_module

def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if parallelism == "OpenMP":
        if solver_type == "static" or solver_type == "Static":
            solver_module_name = "adaptative_remeshing_contact_structural_mechanics_static_solver"
        elif solver_type == "dynamic" or solver_type == "Dynamic":
            solver_module_name = "adaptative_remeshing_contact_structural_mechanics_implicit_dynamic_solver"
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.ContactStructuralMechanicsApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, custom_settings["solver_settings"])

    return solver
