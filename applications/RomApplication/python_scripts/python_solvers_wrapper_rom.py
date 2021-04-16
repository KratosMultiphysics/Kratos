from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = solver_settings["solver_type"].GetString()


    if solver_settings.Has("time_integration_method"):
        time_integration_method = solver_settings["time_integration_method"].GetString()
    else:
        time_integration_method = "implicit" # defaulting to implicit time-integration    

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "transient" or solver_type == "Transient"):
            solver_module_name = "convection_diffusion_transient_rom_solver"
  
        elif (solver_type == "dynamic" or solver_type == "Dynamic"):
            if time_integration_method == "implicit":
                solver_module_name = "structural_mechanics_implicit_dynamic_rom_solver" 
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\""
                raise Exception(err_msg)

        elif solver_type == "static" or solver_type == "Static":
            solver_module_name = "structural_mechanics_static_rom_solver"                       

        elif (solver_type == "stationary" or solver_type == "Stationary"):
            solver_module_name = "convection_diffusion_stationary_rom_solver"
            
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"transient\", \"stationary\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        err_msg =  "The requested parallel type MPI is not yet available!\n"
        raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.RomApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism)

