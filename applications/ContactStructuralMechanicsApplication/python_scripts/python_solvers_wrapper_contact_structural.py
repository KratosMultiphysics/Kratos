from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "dynamic" or solver_type == "Dynamic"):
            solver_module_name = "contact_structural_mechanics_implicit_dynamic_solver"

        elif (solver_type == "static" or solver_type == "Static"):
            solver_module_name = "contact_structural_mechanics_static_solver"
        else:
            raise Exception("The requested solver type is not in the python solvers wrapper")

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        raise Exception("The requested solver type is not in the python solvers wrapper")
        #if (solver_type == "dynamic" or solver_type == "Dynamic"):
            #solver_module_name = "trilinos_contact_structural_mechanics_implicit_dynamic_solver"

        #elif (solver_type == "static" or solver_type == "Static"):
            #solver_module_name = "trilinos_contact_structural_mechanics_static_solver"

        #else:
            #raise Exception("The requested solver type is not in the python solvers wrapper")
    else:
        raise Exception("Parallelism is neither OpenMP nor MPI")

    # Remove settings that are not needed any more
    custom_settings["solver_settings"].RemoveValue("solver_type")
    custom_settings["solver_settings"].RemoveValue("time_integration_method") # does not throw even if the value is not existing

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, custom_settings["solver_settings"])

    return solver
