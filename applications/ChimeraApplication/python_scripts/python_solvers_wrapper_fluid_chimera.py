# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
from importlib import import_module
from KratosMultiphysics.ChimeraApplication import *


def CreateSolverByParameters(model, custom_settings, parallelism):

    # initializing the solver module
    solver_module_name = ""

    solver_type = custom_settings["solver_type"].GetString()

    # Turbulent Chimera Solver
    '''Parallelism check for turbulent flows is not necessary as this is 
       done automatically in RANS Application'''
    if (solver_type == "monolithic_rans_chimera" or solver_type == "MonolithicRANSChimera"):
        # Check for <monolithic> keyword in formulation_name of formulation_settings
        if "monolithic" in custom_settings["formulation_settings"]["formulation_name"].GetString():
            solver_module_name = "rans_solver_vmsmonolithic_chimera"
        else:
            err_msg = "The solver type <{}> is of a monolithic approach.\n".format(solver_type)
            err_msg += "However, the RANS formulation_name is of different approach. "
            err_msg += "Both solver_type and RANS formulation_name must be of monolithic approach."
            raise Exception(err_msg)

    # Laminar Chimera Solver for OpenMP parallelism
    elif (parallelism == "OpenMP"):
        if (solver_type == "monolithic" or solver_type == "Monolithic"):
            solver_module_name = "navier_stokes_solver_vmsmonolithic_chimera"
        elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
            solver_module_name = "navier_stokes_solver_fractionalstep_chimera"

    # Laminar Chimera Solver for MPI parallelism
    elif (parallelism == "MPI"):
        raise Exception("Only OpenMP parallelism supported for <monolithic> and <fractional_step> solver_type")
    else:
        err_msg = "The requested solver type <{}> is not supported. ".format(solver_type)
        err_msg += "Only the following solver types are supported: \n"
        err_msg += "monolithic_rans_chimera\n",
        err_msg += "monolithic\n"
        err_msg += "fractional_step"
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.ChimeraApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, custom_settings)

    return solver


def CreateSolver(model, custom_settings):

    if (not isinstance(model, KratosMultiphysics.Model)):
        raise Exception(
            "input is expected to be provided as a Kratos Model object")

    if (not isinstance(custom_settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
