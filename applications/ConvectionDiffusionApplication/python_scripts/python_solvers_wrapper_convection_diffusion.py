from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = solver_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "transient" or solver_type == "Transient"):
            solver_module_name = "convection_diffusion_transient_solver"

        elif (solver_type == "stationary" or solver_type == "Stationary"):
            solver_module_name = "convection_diffusion_stationary_solver"

        elif (solver_type == "thermally_coupled" or solver_type == "ThermallyCoupled"):
            solver_module_name = "coupled_fluid_thermal_solver"

        elif (solver_type == "thermo_mechanically_coupled" or solver_type == "ThermoMechanicallyCoupled"):
            solver_module_name = "coupled_structural_thermal_solver"

        elif (solver_type == "conjugate_heat_transfer" or solver_type == "ConjugateHeatTransfer"):
            solver_module_name = "conjugate_heat_transfer_solver"

        elif solver_type == "adjoint_stationary":
            solver_module_name = "adjoint_diffusion_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"transient\", \"stationary\", \"thermally_coupled\", \"conjugate_heat_transfer\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        err_msg =  "The requested parallel type MPI is not yet available!\n"
        raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.ConvectionDiffusionApplication.' + solver_module_name
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

