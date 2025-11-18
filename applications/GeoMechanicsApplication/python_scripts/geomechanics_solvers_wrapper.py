import KratosMultiphysics
from importlib import import_module


def CreateSolverByParameters(model, custom_settings, parallelism):
    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type_raw = custom_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        solver_type = solver_type_raw.lower()
        if solver_type in ("u_pw", "geomechanics_u_pw_solver", "twophase"):
            solver_module_name = "geomechanics_U_Pw_solver"

        elif solver_type in ("pw", "geomechanics_pw_solver"):
            solver_module_name = "geomechanics_Pw_solver"

        elif solver_type in ("t", "geomechanics_t_solver"):
            solver_module_name = "geomechanics_T_solver"
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"geomechanics_U_Pw_solver\", \"geomechanics_Pw_solver\", \"geomechanics_T_solver\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\""
        raise Exception(err_msg)

    module_full_name = 'KratosMultiphysics.GeoMechanicsApplication.' + solver_module_name
    return import_module(module_full_name).CreateSolver(model, custom_settings)

def ExtractModelPartNames(process_list, domain_condition_names, root_name, prefix):
    for i in range(process_list.size()):
        process = process_list[i]
        if process.Has("Parameters") and process["Parameters"].Has("model_part_name"):
            model_part_name = process["Parameters"]["model_part_name"].GetString()
            if model_part_name == root_name:
               continue
            if model_part_name.startswith(prefix):
               model_part_name = model_part_name[len(prefix):]
            domain_condition_names.add(model_part_name)

def CreateSolver(model, custom_settings):
    solver_settings = custom_settings["solver_settings"]
    if solver_settings.Has("processes_sub_model_part_list"):
        solver_settings.RemoveValue("processes_sub_model_part_list")
    solver_settings.AddEmptyArray("processes_sub_model_part_list")

    process_lists_to_be_checked = [
        "constraints_process_list",
        "loads_process_list",
        "auxiliary_process_list"
    ]
    domain_condition_names = set()
    root_name = solver_settings["model_part_name"].GetString()
    prefix = root_name + "."
    for process_list_name in process_lists_to_be_checked:
        if custom_settings["processes"].Has(process_list_name):
            ExtractModelPartNames(
               custom_settings["processes"][process_list_name],
               domain_condition_names,
               root_name,
               prefix
        )
    for name in domain_condition_names:
        solver_settings["processes_sub_model_part_list"].Append(KratosMultiphysics.Parameters(f'"{name}"'))

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    return CreateSolverByParameters(model, solver_settings, parallelism)