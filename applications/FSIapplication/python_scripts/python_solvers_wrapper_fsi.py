from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["coupling_solver_settings"]["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["coupling_solver_settings"]["solver_settings"]["solver_type"].GetString()
    coupling_scheme = custom_settings["coupling_solver_settings"]["solver_settings"]["coupling_scheme"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if solver_type == "Partitioned" and coupling_scheme == "DirichletNeumann":
            solver_module_name = "partitioned_fsi_dirichlet_neumann_solver"
        else:
            err_msg =  'The next combination is not available in the FSI Python solvers wrapper.'
            err_msg =  '\t- Requested solver_type: ' + solver_type + '\n'
            err_msg += '\t- Requested coupling_scheme: ' + coupling_scheme + '\n'
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if solver_type == "Partitioned" and coupling_scheme == "DirichletNeumann":
            solver_module_name = "trilinos_partitioned_fsi_dirichlet_neumann_solver"
        else:
            err_msg =  'The next combination is not available in the FSI Python solvers wrapper.'
            err_msg =  '\t- Requested solver_type: ' + solver_type + '\n'
            err_msg += '\t- Requested coupling_scheme: ' + coupling_scheme + '\n'
            raise Exception(err_msg)
    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, custom_settings)

    return solver
