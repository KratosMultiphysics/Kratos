from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_settings.Has("time_integration_method"):
        time_integration_method = solver_settings["time_integration_method"].GetString()
    else:
        time_integration_method = "implicit" # defaulting to implicit time-integration

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "dynamic" or solver_type == "Dynamic"):
            if (time_integration_method == "implicit"):
                solver_module_name = "structural_mechanics_implicit_dynamic_solver"
            elif ( time_integration_method == "explicit"):
                solver_module_name = "structural_mechanics_explicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\", \"explicit\""
                raise Exception(err_msg)

        elif (solver_type == "static" or solver_type == "Static"):
            solver_module_name = "structural_mechanics_static_solver"

        elif (solver_type == "eigen_value"):
            solver_module_name = "structural_mechanics_eigensolver"

        elif (solver_type == "harmonic_analysis"):
            solver_module_name = "structural_mechanics_harmonic_analysis_solver"

        elif (solver_type == "formfinding"):
            solver_module_name = "structural_mechanics_formfinding_solver"

        elif (solver_type == "adjoint_static"):
            solver_module_name = "structural_mechanics_adjoint_static_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\", \"eigen_value\", \"harmonic_analysis\", \"formfinding\", \"adjoint_static\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "dynamic" or solver_type == "Dynamic"):
            if (time_integration_method == "implicit"):
                solver_module_name = "trilinos_structural_mechanics_implicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\""
                raise Exception(err_msg)

        elif (solver_type == "static" or solver_type == "Static"):
            solver_module_name = "trilinos_structural_mechanics_static_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    # Remove settings that are not needed any more
    solver_settings.RemoveValue("solver_type")
    solver_settings.RemoveValue("time_integration_method") # does not throw even if the value is not existing

    module_full = 'KratosMultiphysics.StructuralMechanicsApplication.' + solver_module_name
    solver = __import__(module_full, fromlist=[solver_module_name]).CreateSolver(model, solver_settings)

    return solver


def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
