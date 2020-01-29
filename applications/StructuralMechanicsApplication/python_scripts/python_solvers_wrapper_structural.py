from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

# Other imports
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_settings.Has("time_integration_method"):
        time_integration_method = solver_settings["time_integration_method"].GetString()
    else:
        time_integration_method = "implicit" # defaulting to implicit time-integration

    # Solvers for OpenMP parallelism
    if parallelism == "OpenMP":
        if solver_type == "dynamic" or solver_type == "Dynamic":
            if time_integration_method == "implicit":
                solver_module_name = "structural_mechanics_implicit_dynamic_solver"
            elif time_integration_method == "explicit":
                solver_module_name = "structural_mechanics_explicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\", \"explicit\""
                raise Exception(err_msg)

        elif solver_type == "static" or solver_type == "Static":
            solver_module_name = "structural_mechanics_static_solver"

        elif solver_type == "eigen_value":
            solver_module_name = "structural_mechanics_eigensolver"

        elif solver_type == "harmonic_analysis":
            solver_module_name = "structural_mechanics_harmonic_analysis_solver"

        elif solver_type == "formfinding":
            solver_module_name = "structural_mechanics_formfinding_solver"

        elif solver_type == "adjoint_static":
            solver_module_name = "structural_mechanics_adjoint_static_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\", \"eigen_value\", \"harmonic_analysis\", \"formfinding\", \"adjoint_static\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif parallelism == "MPI":
        if solver_type == "dynamic" or solver_type == "Dynamic":
            if time_integration_method == "implicit":
                solver_module_name = "trilinos_structural_mechanics_implicit_dynamic_solver"
            else:
                err_msg =  "The requested time integration method \"" + time_integration_method + "\" is not in the python solvers wrapper\n"
                err_msg += "Available options are: \"implicit\""
                raise Exception(err_msg)

        elif solver_type == "static" or solver_type == "Static":
            solver_module_name = "trilinos_structural_mechanics_static_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"static\", \"dynamic\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    if solver_settings.Has("contact_settings"):  # This is a contact problem
        kratos_module = "KratosMultiphysics.ContactStructuralMechanicsApplication"
        solver_module_name = "contact_" + solver_module_name

        if solver_settings.Has("use_computing_model_part"):
            if not solver_settings["use_computing_model_part"].GetBool():
                KratosMultiphysics.Logger.PrintWarning('MechanicalSolver', 'Contact Problems must use the computing model part for now, setting it to on')
        else:
            solver_settings.AddEmptyValue("use_computing_model_part")
        solver_settings["use_computing_model_part"].SetBool(True)

    elif solver_settings.Has("mpc_contact_settings"):  # This is a mpc contact problem
        kratos_module = "KratosMultiphysics.ContactStructuralMechanicsApplication"
        solver_module_name = "mpc_contact_" + solver_module_name

        if solver_settings.Has("use_computing_model_part"):
            if not solver_settings["use_computing_model_part"].GetBool():
                KratosMultiphysics.Logger.PrintWarning('MechanicalSolver', 'Contact Problems must use the computing model part for now, setting it to on')
        else:
            solver_settings.AddEmptyValue("use_computing_model_part")
        solver_settings["use_computing_model_part"].SetBool(True)

    else:
        kratos_module = "KratosMultiphysics.StructuralMechanicsApplication"

        if solver_settings.Has("use_computing_model_part"):
            if solver_settings["use_computing_model_part"].GetBool():
                kratos_utils.IssueDeprecationWarning('MechanicalSolver', 'Using ["use_computing_model_part" : True] is deprecated, please remove this setting')
        elif not (solver_settings.Has("problem_domain_sub_model_part_list") or solver_settings.Has("processes_sub_model_part_list")):
            solver_settings.AddEmptyValue("use_computing_model_part").SetBool(False)


    solver = import_module(kratos_module + "." + solver_module_name).CreateSolver(model, solver_settings)

    return solver


def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
