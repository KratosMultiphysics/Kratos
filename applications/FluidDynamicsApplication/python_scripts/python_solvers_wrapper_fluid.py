from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_type == "ale_fluid":
        import navier_stokes_ale_fluid_solver
        return navier_stokes_ale_fluid_solver.CreateSolver(model, solver_settings, parallelism)

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Monolithic"):
            solver_module_name = "navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "navier_stokes_solver_fractionalstep"

        elif (solver_type == "Embedded"):
            solver_module_name = "navier_stokes_embedded_solver"

        elif (solver_type == "Compressible"):
            solver_module_name = "navier_stokes_compressible_solver"

        elif (solver_type == "ConjugateHeatTransfer"):
            solver_module_name = "conjugate_heat_transfer_solver"

        else:
            raise Exception("the requested solver type is not in the python solvers wrapper. Solver type is : " + solver_type)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Monolithic"):
            solver_module_name = "trilinos_navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "trilinos_navier_stokes_solver_fractionalstep"

        elif (solver_type == "Embedded"):
            solver_module_name = "trilinos_navier_stokes_embedded_solver"

        else:
            raise Exception("the requested solver type is not in the python solvers wrapper. Solver type is : " + solver_type)

    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
