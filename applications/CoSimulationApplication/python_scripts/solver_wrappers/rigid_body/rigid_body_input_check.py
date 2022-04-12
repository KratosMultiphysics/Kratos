
import KratosMultiphysics as KM

import json
from importlib import import_module


def _CheckMandatoryInputParameters(project_parameters):

    for key in ["problem_data", "solver_settings", "output_processes", "processes"]:
        if not project_parameters.Has(key):
            msg = 'The key "' + key + '" was not found in the project parameters '
            msg += 'and it is necessary to configure the RigidBodySolver.'
            raise Exception(msg)
    
    for key in ["problem_name", "start_time", "end_time"]:
        if not project_parameters["problem_data"].Has(key):
            msg = '"'+key+'" should be given as par of "problem_data" in the project parameters.'
            raise Exception(msg)

    msg = '"time_step" should be given as par of "time_integration_parameters" '
    msg += 'in "solver_settings" of the project parameters.'
    if not project_parameters["solver_settings"].Has("time_integration_parameters"):
        raise Exception(msg)
    else:
        if not project_parameters["solver_settings"]["time_integration_parameters"].Has("time_step"):
            raise Exception(msg)

    if project_parameters["solver_settings"].Has("model_import_settings"):
        if project_parameters["solver_settings"]["model_import_settings"].Has("input_type"):
            input_type = project_parameters["solver_settings"]["model_import_settings"]["input_type"].GetString()
            if input_type not in ["none","rest"]:
                raise Exception('"input_type" must be either "none" or "rest".')
            if input_type == "rest":
                if not project_parameters["solver_settings"]["model_import_settings"].Has("restart_load_file_label"):
                    raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')


def _ValidateAndAssignRigidBodySolverDefaults(solver_settings):
    
    default_solver_settings = KM.Parameters('''{
        "domain_size":3,
        "echo_level":0,
        "buffer_size":3,
        "model_import_settings":{
            "input_type":"none",
            "input_filename":"Main",
            "restart_load_file_label":"0.0",
            "load_restart_files_from_folder":true,
            "input_output_path":"restart"
        },
        "time_integration_parameters":{
            "rho_inf": 0.16,
            "time_step": 0.1
        },
        "active_dofs":[]
    }''')

    solver_settings.RecursivelyValidateAndAssignDefaults(default_solver_settings)

    if solver_settings["domain_size"].GetInt() not in [2, 3]:
        raise Exception("The domain size can only be 2 or 3.")
    
    if solver_settings["buffer_size"].GetInt() < 1:
        raise Exception("The buffer size needs to be equal or bigger than 1.")

    return solver_settings


def _ValidateAndAssignDofDefaults(dof_settings, available_dofs):

    if dof_settings.size() == 0:
        msg = 'At least an active degree of freedom is needed to use the solver '
        msg += ' and none where provided in "active_dofs".'
        raise Exception(msg)

    default_dof_settings = KM.Parameters('''{
        "dof":"needs_to_be_given",
        "constrained": false,
        "system_parameters":{
            "mass"      : 1.0,
            "stiffness" : 1.0,
            "damping"   : 0.0,
            "modulus_self_weight": 0.0
        }
    }''')

    dof_settings_processed = {}
    active_dofs = []

    # Read entries from the project parameters and check are not repeated
    for entry in dof_settings:
        dof_data = entry.Clone()
        dof = dof_data["dof"].GetString()
        if dof not in available_dofs:
            msg = 'The degree of freedom "'+dof+'" is not among the available ones. '
            msg += 'Select one of the following: '+str(available_dofs)[1,-1]
            raise Exception(msg)
        if dof in dof_settings_processed:
            msg = 'The degree of freedom "'+dof+'" was repeated in the list "active dofs".'
            raise Exception(msg)
        else:
            dof_data.RecursivelyValidateAndAssignDefaults(default_dof_settings)
            dof_settings_processed[dof] = dof_data

    # Fill dofs which were not mentioned explicitly with the default values
    for dof in available_dofs:
        if dof not in dof_settings_processed:
            dof_settings_processed[dof] = default_dof_settings
        else:
            active_dofs.append(dof)

    return dof_settings_processed, active_dofs


def _CreateListOfProcesses(model, parameters, main_model_part):

    # Create all the processes stated in the project parameters
    # TODO: No need to convert it to a json to manipulate it
    if main_model_part.ProcessInfo[KM.IS_RESTARTED]:
        process_types = ["gravity", "boundary_conditions_process_list", "auxiliar_process_list"]
    else:
        process_types = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
    parameters_json = json.loads(parameters.WriteJsonString())
    list_of_processes = []
    # TODO: Is this usually a mandatory input?
    if "processes" in parameters_json:
        for process_type in process_types:
            if process_type in parameters_json["processes"]:
                for process in parameters_json["processes"][process_type]:
                    python_module = process["python_module"]
                    kratos_module = process["kratos_module"]
                    process_module = import_module(kratos_module + "." + python_module)
                    process_settings = KM.Parameters(json.dumps(process))
                    list_of_processes.append(process_module.Factory(process_settings, model))
    
    return list_of_processes


def _CreateListOfOutputProcesses(model, parameters):

    # Create all the processes stated in the project parameters
    # TODO: No need to convert it to a json to manipulate it
    parameters_json = json.loads(parameters.WriteJsonString())
    list_of_output_processes = []
    # TODO: Is this usually a mandatory input?
    if "output_processes" in parameters_json:
        for process in parameters_json["output_processes"]:
            python_module = process["python_module"]
            kratos_module = process["kratos_module"]
            process_module = import_module(kratos_module + "." + python_module)
            process_settings = KM.Parameters(json.dumps(process))
            list_of_output_processes.append(process_module.Factory(process_settings, model))
    
    return list_of_output_processes