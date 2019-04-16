from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolverByParameters(model, solver_settings, parallelism):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = solver_settings["solver_type"].GetString()
    coupling_scheme = solver_settings["coupling_scheme"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Partitioned" or solver_type == "partitioned"):
            if (coupling_scheme == "DirichletNeumann" or coupling_scheme == "dirichlet_neumann"):
                solver_module_name = "partitioned_fsi_dirichlet_neumann_solver"
            else:
                err_msg = 'Requested coupling_scheme: ' + coupling_scheme + ' is not available.'
                raise Exception(err_msg)
        else:
            err_msg = 'Requested solver_type: ' + solver_type + ' is not available.'
            raise Exception(err_msg)
    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Partitioned" or solver_type == "partitioned"):
            if (coupling_scheme == "DirichletNeumann" or coupling_scheme == "dirichlet_neumann"):
                solver_module_name = "trilinos_partitioned_fsi_dirichlet_neumann_solver"
            else:
                err_msg = 'Requested coupling_scheme: ' + coupling_scheme + ' is not available.'
                raise Exception(err_msg)
        else:
            err_msg = 'Requested solver_type: ' + solver_type + ' is not available.'
            raise Exception(err_msg)
    else:
        err_msg = "Parallelism is neither OpenMP nor MPI."
        raise Exception(err_msg)

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism)