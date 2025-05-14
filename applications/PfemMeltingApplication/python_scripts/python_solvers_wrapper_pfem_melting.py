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
        # Coupled CFD-thermal solvers (volume coupling by Boussinesq approximation)
        if (solver_type == "thermally_coupledpfem2" or solver_type == "ThermallyCoupledPfem2"):

            solver_module_name = "coupled_fluid_thermal_pfem2solver"

        elif (solver_type == "thermally_coupled_no_remesh"):

            solver_module_name = "coupled_fluid_thermal_solverwithoutmeshgeneration"

        elif (solver_type == "WithoutmeshSolid"):

            solver_module_name = "coupled_solid_thermal_solverwithoutmeshgeneration"
            
            
        elif (solver_type == "WithmeshSolid"):
            
            #solver_module_name = "coupled_solid_thermal_solverwithmeshgeneration"
            
            #solver_module_name = "coupled_fsi_solver" ###fsi limpio
            
            solver_module_name = "coupled_cavity_solver" ###para eugenio

            
            #solver_module_name = "coupled_fsi_solverwithmeshgeneration_complete" ###con este trabaje el ejemplo debug

        # Coupled mechanical-thermal solver
        # Wrong solver check
        else:
            err_msg =  "The requested solver type {} is not in the python solvers wrapper\n".format(solver_type)
            err_msg += "Available options are: \"transient\", \"stationary\", \"thermally_coupled\", \"thermo_mechanically_coupled\", \"conjugate_heat_transfer\" and \"adjoint_stationary\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        err_msg =  "The requested parallel type MPI is not yet available!\n"
        raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.PfemMeltingApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")
    lllllllllllllllllllllllllllll
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism)
