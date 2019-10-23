from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from importlib import import_module

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Monolithic"):
            solver_module_name = "adjoint_vmsmonolithic_solver"
        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Monolithic"):
            solver_module_name = "trilinos_adjoint_vmsmonolithic_solver"
        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")
    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    module_full = 'KratosMultiphysics.FluidDynamicsApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, custom_settings["solver_settings"])

    return solver
