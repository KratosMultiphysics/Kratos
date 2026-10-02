
import KratosMultiphysics as KM

import json
from importlib import import_module


"""
This file checks the input for the RigidBodySolver and raises exceptions if necessary.
"""


def _CheckMandatoryInputParameters(project_parameters):
    # Function to raise exceptions if some mandatory parameters were not given

    # The project parameters need to have at least these base fields
    for key in ["problem_data", "solver_settings", "output_processes", "processes"]:
        if not project_parameters.Has(key):
            msg = f'The key "{key}" was not found in the project parameters '
            msg += 'and it is necessary to configure the RigidBodySolver.'
            raise Exception(msg)
    
    # Some subfields of "problem_data" must be included
    for key in ["problem_name", "start_time", "end_time"]:
        if not project_parameters["problem_data"].Has(key):
            msg = '"'+key+'" should be given as par of "problem_data" in the project parameters.'
            raise Exception(msg)

    # Check that the time step is given
    msg = '"time_step" should be given as par of "time_integration_parameters" '
    msg += 'in "solver_settings" of the project parameters.'
    if not project_parameters["solver_settings"].Has("time_integration_parameters"):
        raise Exception(msg)
    else:
        if not project_parameters["solver_settings"]["time_integration_parameters"].Has("time_step"):
            raise Exception(msg)

    # Check that "input_type", if needed, is among the available ones
    if project_parameters["solver_settings"].Has("model_import_settings"):
        if project_parameters["solver_settings"]["model_import_settings"].Has("input_type"):
            input_type = project_parameters["solver_settings"]["model_import_settings"]["input_type"].GetString()
            if input_type not in ["none","rest"]:
                raise Exception('"input_type" must be either "none" or "rest".')
            if input_type == "rest":
                if not project_parameters["solver_settings"]["model_import_settings"].Has("restart_load_file_label"):
                    raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')


def _ValidateAndAssignRigidBodySolverDefaults(solver_settings):
    # Function to check and fill the solver settings
    
    # Omitted fields will be filled with these default parameters
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

    # The domain size needs to be either 2 or 3
    if solver_settings["domain_size"].GetInt() not in [2, 3]:
        raise Exception("The domain size can only be 2 or 3.")
    if solver_settings["domain_size"].GetInt() == 2:
        msg = 'The 2D version of the solver is yet to be implemented. Use 3 as '
        msg += '"domain_size" and activate only the necessary degrees of freedom instead.'
        raise Exception(msg)
    
    # The generalized alpha time integration scheme needs at least a buffer size of one
    if solver_settings["buffer_size"].GetInt() < 1:
        raise Exception("The buffer size needs to be equal or bigger than 1.")

    return solver_settings


def _ValidateAndAssignDofDefaults(dof_settings, available_dofs):
    # Function to check and fill the structural parameters of each degree of freedom (DOF)

    # Check that at least one DOF was mentioned. If not, tthere is nothing to be simulated
    if dof_settings.size() == 0:
        msg = 'At least an active degree of freedom is needed to use the solver '
        msg += ' and none where provided in "active_dofs".'
        raise Exception(msg)

    # Omitted structural parameters will be filled with these for each DOF
    default_dof_settings = KM.Parameters('''{
        "dof":"needs_to_be_given",
        "constrained": false,
        "system_parameters":{
            "mass"      : 1.0,
            "stiffness" : 1.0,
            "damping"   : 0.0
        }
    }''')

    # Initialize variables to be filled and returned
    dof_settings_processed = {}
    active_dofs = []

    # Read entries from the project parameters and check they are not repeated
    for entry in dof_settings.values():
        # Clone the data so the modifications don't affect previously processed entries
        dof_data = entry.Clone()
        # Check that the selected DOF is among the available ones
        dof = dof_data["dof"].GetString()
        if dof not in available_dofs:
            msg = 'The degree of freedom "'+dof+'" is not among the available ones. '
            msg += 'Select one of the following: '+str(available_dofs)[1:-1]
            raise Exception(msg)
        # Repeated DOFs are not allowed since they can create contradictions
        if dof in dof_settings_processed:
            msg = 'The degree of freedom "'+dof+'" was repeated in the list "active dofs".'
            raise Exception(msg)
        else:
            # Fill with defaults and save in the variable to be outputted
            dof_data.RecursivelyValidateAndAssignDefaults(default_dof_settings)
            dof_settings_processed[dof] = dof_data

    # Fill dofs which were not mentioned explicitly with the default values
    for dof in available_dofs:
        if dof not in dof_settings_processed:
            dof_settings_processed[dof] = default_dof_settings
        else:
            active_dofs.append(dof)

    # Return the data needed by the solver
    return dof_settings_processed, active_dofs


def _CreateListOfProcesses(model, parameters, main_model_part):
    # This function creates all the processes stated in the project parameters

    # If it is a restart, ignore the "initial_conditions_process_list", since they are given by the restart file
    if main_model_part.ProcessInfo[KM.IS_RESTARTED]:
        process_types = ["gravity", "boundary_conditions_process_list", "auxiliar_process_list"]
    else:
        process_types = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]

    # Import the processes and save them in a list
    list_of_processes = []
    if parameters.Has("processes"):
        for process_type in process_types:
            if parameters["processes"].Has(process_type):
                for process in parameters["processes"][process_type].values():
                    if process.Has("name"):
                        registry_entry = process["name"].GetString()
                        if KM.Registry.HasItem(registry_entry):
                            # Get already stored prototype
                            if KM.Registry.HasItem(f"{registry_entry}.Prototype"):
                                prototype = KM.Registry[f"{registry_entry}.Prototype"]
                                instance = prototype.Create(model, process["Parameters"])
                            # Get prototype from stored Python module
                            elif KM.Registry.HasItem(f"{registry_entry}.ModuleName"):
                                class_name = registry_entry.split(".")[-1]
                                module_name = KM.Registry[f"{registry_entry}.ModuleName"]
                                module = import_module(module_name)
                                if hasattr(module, class_name):
                                    prototype = getattr(module, class_name)
                                    instance = prototype(model, process["Parameters"])
                                else:
                                    #TODO: In here we're assuming that the registry last key is the class name
                                    #TODO: We should enforce this. Now an error happens but as we populate we should throw a warning and search for a ClassName item
                                    err_msg = f"The '{class_name}' class name cannot be found within the '{module_name}' module."
                                    raise Exception(err_msg)
                            else:
                                err_msg = f"Registry process '{registry_entry}' cannot be constructed."
                                raise Exception(err_msg)
                        else:
                            KM.Logger.PrintWarning(f"Asking to construct the non-registered 'name' '{registry_entry}'.")
                    else:
                        python_module = process["python_module"].GetString()
                        kratos_module = process["kratos_module"].GetString()
                        process_module = import_module(kratos_module + "." + python_module)
                        instance = process_module.Factory(process, model)
                    # Construct the process from the obtained prototype and append it to the list
                    list_of_processes.append(instance)
    
    return list_of_processes

def _CreateListOfOutputProcesses(model, parameters):
    # This function creates all the output processes stated in the project parameters

    # Import the processes and save them in a list
    list_of_output_processes = []
    if parameters.Has("output_processes"):
        for process in parameters["output_processes"].values():
            python_module = process["python_module"].GetString()
            kratos_module = process["kratos_module"].GetString()
            process_module = import_module(kratos_module + "." + python_module)
            list_of_output_processes.append(process_module.Factory(process, model))
    
    return list_of_output_processes