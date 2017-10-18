from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CheckForParallelism(settings):
    if not settings["problem_data"].Has("parallel_type"):
        raise Exception("\"parallel_type\" is not specified in the settings, please provide it")

def CheckForSolverType(settings):
    if not settings["solver_settings"].Has("solver_type"):
        raise Exception("\"solver_type\" is not specified in the settings, please provide it")

def CheckForAnalysisType(settings):
    if not settings["solver_settings"].Has("analysis_type"):
        raise Exception("\"analysis_type\" is not specified in the settings, please provide it")

def CheckForTimeIntegrationMethod(settings):
    if not settings["solver_settings"].Has("time_integration_method"):
        raise Exception("\"time_integration_method\" is not specified in the settings, please provide it")
    

def CreateSolver(main_model_part, custom_settings):

    if (type(main_model_part) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    CheckForParallelism(custom_settings)
    CheckForSolverType(custom_settings)

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Dynamic"):
            CheckForAnalysisType(custom_settings)
            CheckForTimeIntegrationMethod(custom_settings)
            time_integration_method = custom_settings["solver_settings"]["time_integration_method"].GetString()
            if (time_integration_method == "Implicit"):
                solver_module_name = "structural_mechanics_implicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"Implicit\""
                raise Exception(err_msg)

        elif (solver_type == "Static"):
            CheckForAnalysisType(custom_settings)
            solver_module_name = "structural_mechanics_static_solver"

        elif (solver_type == "EigenValue"):
            solver_module_name = "structural_mechanics_eigensolver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"Static\", \"Dynamic\", \"EigenValue\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Dynamic"):
            CheckForAnalysisType(custom_settings)
            CheckForTimeIntegrationMethod(custom_settings)
            time_integration_method = custom_settings["solver_settings"]["time_integration_method"].GetString()
            if (time_integration_method == "Implicit"):
                solver_module_name = "trilinos_structural_mechanics_implicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"Implicit\""
                raise Exception(err_msg)

        elif (solver_type == "Static"):
            CheckForAnalysisType(custom_settings)
            solver_module_name = "trilinos_structural_mechanics_static_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"Static\", \"Dynamic\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallel_type + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(main_model_part, custom_settings["solver_settings"])

    return solver
