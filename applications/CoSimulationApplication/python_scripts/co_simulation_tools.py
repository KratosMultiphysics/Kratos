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
