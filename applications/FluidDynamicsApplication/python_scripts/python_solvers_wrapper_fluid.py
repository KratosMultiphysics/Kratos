from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def ConstructPythonSolver(main_model_part, ProjectParameters):

    if (type(main_model_part) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if (type(ProjectParameters) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = ProjectParameters["problem_data"]["parallel_type"].GetString()
    solver_type = ProjectParameters["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Monolithic"):
            solver_module_name = "navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "navier_stokes_solver_fractionalstep"

        elif (solver_type == "Embedded"):
            solver_module_name = "navier_stokes_embedded_solver"

        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Monolithic"):
            solver_module_name = "trilinos_navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "trilinos_navier_stokes_solver_fractionalstep"

        elif (solver_type == "Embedded"):
            solver_module_name = "trilinos_navier_stokes_embedded_solver"

        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")
    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

    return solver
