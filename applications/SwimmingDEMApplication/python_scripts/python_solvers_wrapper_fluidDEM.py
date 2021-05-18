from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_type == "ale_fluid":
        import navier_stokes_ale_fluid_solver
        return navier_stokes_ale_fluid_solver.CreateSolver(model, solver_settings, parallelism)

    # Solvers for OpenMP or MPI parallelism
    if (parallelism == "OpenMP"):
        solver_module_name = "navier_stokes_solver_vms_monolithic_DEMCoupled"
    else:
        raise Exception("parallelism is not OpenMP")

    module_full = 'KratosMultiphysics.SwimmingDEMApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)