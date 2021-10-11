
import KratosMultiphysics

import json

# TODO: Check this function
# Check also that the input value is not higher than the buffer size
def _CheckBufferId(buffer_idx, identifier):
    if identifier not in ["DISPLACEMENT", "VELOCITY","ACCELERATION"] and buffer_idx != 0:
        msg = 'The buffer_idx can only be 0 for the variable "' + identifier + '".'
        raise Exception(msg)

def _CheckActiveDofs(parameters, available_dofs):

    active_dofs = []
    for dof in parameters["active_dofs"]:

        if dof not in available_dofs:
            msg = 'The degree of freedom "' + dof + '" is not among the available ones. '
            msg += 'Chose one of the following: ' + str(available_dofs)[1:-1] + '.'
            raise Exception(msg)

        if dof in active_dofs:
            msg = 'The active degree of freedom "' + dof + '" is repeated in the project parameters. '
            raise Exception(msg)
        active_dofs.append(dof)
    
    if len(active_dofs) == 0:
        msg = 'At least an active degree of freedom is needed to use the solver '
        msg += ' and none where provided in "active_dofs".'
        raise Exception(msg)
    
    return active_dofs

def _CheckMandatoryInputParameters(parameters):

    for key in ["active_dofs", "solution_parameters"]:
        if key not in parameters:
            msg = 'The key "' + key + '" was not found in the project parameters '
            msg += 'and it is necessary to configure the RigidBodySolver.'
            raise Exception(msg)
    
    msg = '"time_step" should be given as par of "time_integration_parameters" '
    msg += 'in "solution_parameters".'
    if "time_integration_parameters" not in parameters["solution_parameters"]:
        raise Exception(msg)
    else:
        if "time_step" not in parameters["solution_parameters"]["time_integration_parameters"]:
            raise Exception(msg)

def _ValidateAndAssignRigidBodySolverDefaults(parameters, available_dofs):

    default_single_dof_parameters = KratosMultiphysics.Parameters('''{
        "constrained": false,
        "system_parameters":{
            "mass"      : 1.0,
            "stiffness" : 1.0,
            "damping"   : 0.0,
            "modulus_self_weight": 0.0
        },
        "initial_conditions":{
            "displacement"  : 0.0,
            "velocity"      : 0.0,
            "load_impulse"  : 0.0
        }
    }''')
    default_solution_parameters = KratosMultiphysics.Parameters('''{
        "time_integration_parameters":{
            "rho_inf"   : 0.16,
            "start_time": 0.0,
            "time_step" : 0.05
        },
        "solver_parameters":{
            "buffer_size"   : 2
        },
        "output_parameters":{
            "write_output_files": true,
            "file_path" : "results/rigid_body"
        },
        "restart_parameters":{
            "load_restart_file": false,
            "input_path":   "mdof_solver",
            "input_file_name": "",
            "write_restart_file": false,
            "restart_save_frequency": 0.0,
            "restart_control_type": "time",
            "save_restart_files_in_folder"   : true,
            "output_path"                  : "mdof_solver",
            "max_files_to_keep"            : -1
                }
    }''')

    dof_parameters = {}
    for dof in available_dofs:

        if dof not in parameters["active_dofs"]:
            parameters["active_dofs"][dof] = {}

        single_dof_parameters = parameters["active_dofs"][dof]

        single_dof_parameters = KratosMultiphysics.Parameters(json.dumps(single_dof_parameters))
        single_dof_parameters.RecursivelyValidateAndAssignDefaults(default_single_dof_parameters)
        
        dof_parameters[dof] = single_dof_parameters
    
    solution_parameters = KratosMultiphysics.Parameters(json.dumps(parameters["solution_parameters"]))
    solution_parameters.RecursivelyValidateAndAssignDefaults(default_solution_parameters)

    return dof_parameters, solution_parameters