
import KratosMultiphysics as KM

import json


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
        "output_parameters":{
            "write_output_files": false,
            "file_path" : "results"
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