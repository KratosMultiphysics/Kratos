# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# Importing Kratos
import KratosMultiphysics

# Other imports
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    # Solvers available
    if solver_type == "pfem_fluid_solver" or solver_type == "PfemFluid":
        solver_module_name = "pfem_fluid_solver"

    elif solver_type == "pfem_fluid_nodal_integration_solver" or solver_type == "PfemFluidNodalIntegration":
        solver_module_name = "pfem_fluid_nodal_integration_solver"

    else:
        err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
        err_msg += "Available options are: \"pfem_fluid_solver\", \"pfem_fluid_nodal_integration_solver\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.PfemFluidDynamicsApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver


def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
