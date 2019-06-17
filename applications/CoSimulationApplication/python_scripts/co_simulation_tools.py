import json
cs_data_structure = None

## ImportDataStructure : Imports the data structure which is specified in the parameters file
#
#  @param parameters_file_name   The JSON file name which contains the settings for the co-simulation
def ImportDataStructure(parameters_file_name):
    global cs_data_structure
    if cs_data_structure is None:
        import json
        with open(parameters_file_name,'r') as parameter_file:
            parameters = json.load(parameter_file) # This is for getting the flag for database
            data_structure_name = 'KratosMultiphysics' # Kratos is default
            if 'data_structure' in parameters['problem_data'].keys():
                data_structure_name = parameters['problem_data']['data_structure']
                if not data_structure_name in ['KratosMultiphysics', 'pyKratos']:
                    raise Exception('"data_structure" needs to be "KratosMultiphysics" or "pyKratos"')

            # Initialize cs_data_structure and import corresponding module
            cs_data_structure = __import__(data_structure_name)

    return cs_data_structure


def cs_print_info(label, *args):
    cs_data_structure.Logger.PrintInfo(label, " ".join(map(str,args)))

def cs_print_warning(*args):
    cs_data_structure.Logger.PrintWarning(label, " ".join(map(str,args)))


def CreatePredictors(predictor_settings_list, solvers, parent_solver_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.predictor_factory import CreatePredictor
    predictors = []
    for predictor_settings in predictor_settings_list:
        solver = solvers[predictor_settings["solver"].GetString()]
        if not predictor_settings.Has("echo_level"):
            predictor_settings.AddEmptyValue("echo_level").SetInt(parent_solver_echo_level)
        predictors.append(CreatePredictor(predictor_settings, solver))
    return predictors

def CreateConvergenceAccelerators(convergence_accelerator_settings_list, solvers, parent_solver_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator
    convergence_accelerators = []
    for conv_acc_setting in convergence_accelerator_settings_list:
        solver = solvers[conv_acc_setting["solver"].GetString()]
        if not conv_acc_setting.Has("echo_level"):
            conv_acc_setting.AddEmptyValue("echo_level").SetInt(parent_solver_echo_level)
        convergence_accelerators.append(CreateConvergenceAccelerator(conv_acc_setting, solver))

    return convergence_accelerators

def CreateConvergenceCriteria(convergence_criteria_settings_list, solvers, parent_solver_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.convergence_criteria_factory import CreateConvergenceCriteria
    convergence_criteria = []
    for conv_crit_setting in convergence_criteria_settings_list:
        solver = solvers[conv_crit_setting["solver"].GetString()]
        if not conv_crit_setting.Has("echo_level"):
            conv_crit_setting.AddEmptyValue("echo_level").SetInt(parent_solver_echo_level)
        convergence_criteria.append(CreateConvergenceCriteria(conv_crit_setting, solver))

    return convergence_criteria

def CreateCouplingOperations(coupling_operations_settings_dict, solvers, parent_solver_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.coupling_operation_factory import CreateCouplingOperation
    coupling_operations = {}
    for coupling_operation_name, coupling_operation_setting in coupling_operations_settings_dict.items():
        coupling_operations[coupling_operation_name] = CreateCouplingOperation(coupling_operation_setting, solvers)
        if not coupling_operation_setting.Has("echo_level"):
            coupling_operations[-1].SetEchoLevel(parent_solver_echo_level)

    return coupling_operations

def CreateDataTransferOperators(data_transfer_operators_settings_dict, parent_solver_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.data_transfer_operator_factory import CreateDataTransferOperator
    data_transfer_operators = {}
    for data_transfer_operators_name, data_transfer_operators_setting in data_transfer_operators_settings_dict.items():
        data_transfer_operators[data_transfer_operators_name] = CreateDataTransferOperator(data_transfer_operators_setting)
        if not data_transfer_operators_setting.Has("echo_level"):
            data_transfer_operators[-1].SetEchoLevel(parent_solver_echo_level)

    return data_transfer_operators


def SettingsTypeCheck(settings):
    if not isinstance(settings, cs_data_structure.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
















# def CheckCoSimulationSettingsAndAssignDefaults(co_simulation_settings):
#     # TODO check if the data is consitently defined! => maybe do at another place though...
#     # - input in one is output in another
#     # - one IO is defined for each data_name
#     # - if the same data is defined multiple times
#     # - check if data format has been specified

#     def CheckIfEntryExists(key, dict_name, dict_to_check):
#         if not dict_to_check.Has(key):
#             raise Exception('"' + key + '" not existing in "' + dict_name + '"!')

#     def CheckSolverSetting(solvers):
#         solver_data_dict = {}
#         data_default = cs_data_structure.Parameters(r'''{
#             "data_identifier" : "UNSPECIFIED",
#             "data_format"     : "UNSPECIFIED",
#             "geometry_name"   : "undefined",
#             "type_of_quantity": "undefined"
#         }''')
#         for solver_name, solver_definition in solvers.items():
#             CheckIfEntryExists("data", solver_name, solver_definition)
#             solver_data_dict[solver_name] = solver_definition["data"].keys()
#             for data_name, data_def in solver_definition["data"].items():
#                 data_def.ValidateAndAssignDefaults(data_default)
#                 keys_to_pop = []
#                 for data_key, data_val in data_def.items():
#                     if data_val == "UNSPECIFIED":
#                         err_msg  = '"' + data_key + '" not specified for data "'
#                         err_msg += data_name + '" of solver "' + solver_name + '"'
#                         raise Exception(err_msg)
#                     elif data_val == "undefined":
#                         keys_to_pop.append(data_key)
#                 for key in keys_to_pop: # removing the undefined values again
#                     data_def.pop(key)

#         return solver_data_dict

#     def CheckCouplingLoop(coupling_sequence, solver_data_dict):
#         def CheckInputDataList(input_data_list, solver_name, solver_data_dict):
#             input_data_list_default = {
#                 "from_solver" : "UNSPECIFIED",
#                 "data_name"   : "UNSPECIFIED",
#                 "io_settings" : { }
#             }
#             for data in input_data_list:
#                 ValidateAndAssignDefaults(input_data_list_default, data)
#                 from_solver = data["from_solver"]
#                 data_name   = data["data_name"]
#                 if from_solver == "UNSPECIFIED":
#                     err_msg  = '"from_solver" not specified for solver "'
#                     err_msg += solver_name + '" in "coupling_sequence"!'
#                     raise Exception(err_msg)
#                 if data_name == "UNSPECIFIED":
#                     err_msg  = '"data_name" not specified for solver "'
#                     err_msg += solver_name + '" in "coupling_sequence"!'
#                     raise Exception(err_msg)
#                 if from_solver not in solver_data_dict:
#                     raise Exception('the solver with name "' + from_solver + '" was not specified under "solvers"!')
#                 if data_name not in solver_data_dict[solver_name]: # checking myself
#                     err_msg  = 'the solver with name "' + solver_name
#                     err_msg += '" does not have the data "' + data_name + '"!'
#                     raise Exception(err_msg)
#                 if data_name not in solver_data_dict[from_solver]: # checking the other solver
#                     err_msg  = 'the solver with name "' + from_solver
#                     err_msg += '" does not have the data "' + data_name + '"!'
#                     raise Exception(err_msg)

#         def CheckOutputDataList(output_data_list, solver_name, solver_data_dict):
#             output_data_list_default = {
#                 "to_solver" : "UNSPECIFIED",
#                 "data_name"   : "UNSPECIFIED",
#                 "io_settings" : { }
#             }
#             for data in output_data_list:
#                 ValidateAndAssignDefaults(output_data_list_default, data)
#                 to_solver = data["to_solver"]
#                 data_name = data["data_name"]
#                 if to_solver == "UNSPECIFIED":
#                     err_msg  = '"to_solver" not specified for solver "'
#                     err_msg += solver_name + '" in "coupling_sequence"!'
#                     raise Exception(err_msg)
#                 if data_name == "UNSPECIFIED":
#                     err_msg  = '"data_name" not specified for solver "'
#                     err_msg += solver_name + '" in "coupling_sequence"!'
#                     raise Exception(err_msg)
#                 if to_solver not in solver_data_dict:
#                     raise Exception('the solver with name "' + to_solver + '" was not specified under "solvers"!')
#                 if data_name not in solver_data_dict[solver_name]: # checking myself
#                     err_msg  = 'the solver with name "' + solver_name
#                     err_msg += '" does not have the data "' + data_name + '"!'
#                     raise Exception(err_msg)
#                 if data_name not in solver_data_dict[to_solver]: # checking the other solver
#                     err_msg  = 'the solver with name "' + to_solver
#                     err_msg += '" does not have the data "' + data_name + '"!'
#                     raise Exception(err_msg)

#         solver_default = {
#             "name" : "UNSPECIFIED",
#             "input_data_list"  : [],
#             "output_data_list" : [],
#             "input_coupling_start_time" : 0.0,
#             "output_coupling_start_time" : 0.0
#         }
#         for solver_i, solver in enumerate(coupling_sequence):
#             ValidateAndAssignDefaults(solver_default, solver)
#             solver_name = solver["name"]
#             if solver_name == "UNSPECIFIED":
#                 raise Exception('"name" not specified for solver ' + str(solver_i) + ' in "coupling_sequence"!')
#             if solver_name not in solver_data_dict:
#                 raise Exception('the solver with name "' + solver_name + '" was not specified under "solvers"!')
#             CheckInputDataList(solver["input_data_list"], solver_name, solver_data_dict)
#             CheckOutputDataList(solver["output_data_list"], solver_name, solver_data_dict)

#     problem_data_defaults = cs_data_structure.Parameters(r'''{
#         "start_time"    : -1.0,
#         "end_time"      : -1.0,
#         "echo_level"    : 0,
#         "parallel_type" : "OpenMP",
#         "print_colors"  : false,
#         "flush_stdout"  : false
#     }''')

#     # checking "problem_data"
#     CheckIfEntryExists("problem_data", "co_simulation_settings", co_simulation_settings)
#     problem_data = co_simulation_settings["problem_data"]

#     co_simulation_settings["problem_data"].ValidateAndAssignDefaults(problem_data_defaults)

#     start_time = problem_data["start_time"].GetDouble()
#     if start_time < 0.0:
#         raise Exception('"start_time" has to be >= 0.0!')

#     end_time = problem_data["end_time"].GetDouble()
#     if end_time < 0.0:
#         raise Exception('"end_time" has to be >= 0.0!')
#     if start_time > end_time:
#         raise Exception('"end_time" has to be > "start_time"!')

#     # checking "solver_settings"
#     CheckIfEntryExists("solver_settings", "co_simulation_settings", co_simulation_settings)
#     solver_settings = co_simulation_settings["solver_settings"]

#     CheckIfEntryExists("solver_type", "solver_settings", solver_settings)

#     # This is only existing in a coupled simulation
#     if solver_settings.Has("coupling_sequence"):
#         solver_data_dict = CheckSolverSetting(solver_settings["solvers"])
#         CheckIfEntryExists("solvers", "solver_settings", solver_settings)
#         CheckCouplingLoop(solver_settings["coupling_sequence"], solver_data_dict)



# def ValidateAndAssignDefaults(defaults, settings, recursive=False):
#     for key, val in settings.items():
#         # check if the current entry also exists in the defaults
#         if not key in defaults.keys():
#             err_msg  = 'The item with name "' + key + '" is present in this '
#             err_msg += 'settings\nbut NOT in the defaults!\n'
#             err_msg += 'settings are:\n'
#             err_msg += json.dumps(settings, indent=4)
#             err_msg += '\ndefaults are:\n'
#             err_msg += json.dumps(defaults, indent=4)
#             raise Exception(err_msg)

#         # check if the type is the same in the defaults
#         if type(settings[key]) != type(defaults[key]):
#             err_msg  = 'The type of the item with name "' + key + '" (type: "'
#             err_msg += str(type(settings[key]).__name__)+'") in this '
#             err_msg += 'settings\nis NOT the same as in the defaults (type: "'
#             err_msg += str(type(defaults[key]).__name__)+'")!\n'
#             err_msg += 'settings are:\n'
#             err_msg += json.dumps(settings, indent=4)
#             err_msg += '\ndefaults are:\n'
#             err_msg += json.dumps(defaults, indent=4)
#             raise Exception(err_msg)

#     # loop the defaults and add the missing entries
#     for key_d, val_d in defaults.items():
#         if key_d not in settings: # add the default in case the setting is not present
#             settings[key_d] = val_d
#         elif recursive and type(val_d) is dict:
#             RecursivelyValidateAndAssignDefaults(val_d, settings[key_d])

# def RecursivelyValidateAndAssignDefaults(defaults, settings):
#     ValidateAndAssignDefaults(defaults, settings, recursive=True)

# def csprint(*args):
#     cs_data_structure.Logger.PrintInfo(" ".join(map(str,args)))

# def solverprint(solver_name, *args):
#     csprint(yellow(solver_name + ":"), *args)

# def couplingsolverprint(solver_name, *args):
#     csprint(darkyellow(solver_name + ":"), *args)

# def classprint(solver_name, *args):
#     csprint(magenta(solver_name + ":"), *args)