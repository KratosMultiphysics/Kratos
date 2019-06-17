from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

"""This file contains some functionalities that are commonly used in CoSimulation"""

import json
cs_data_structure = None # global variable that contains the data-structure (Kratos or pyKratos)

def ImportDataStructure(parameters_file_name):
    """Imports the data structure which is specified in the parameters file"""

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

def cs_print_warning(label, *args):
    cs_data_structure.Logger.PrintWarning(label, " ".join(map(str,args)))


def SettingsTypeCheck(settings):
    if not isinstance(settings, cs_data_structure.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")


def AddEchoLevelToSettings(settings, echo_level):
    echo_level_params = cs_data_structure.Parameters("""{
        "echo_level" : """ + str(echo_level) + """
    }""")
    settings.AddMissingParameters(echo_level_params)


def CreatePredictors(predictor_settings_list, solvers, parent_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.predictor_factory import CreatePredictor
    predictors = []
    for predictor_settings in predictor_settings_list:
        solver = solvers[predictor_settings["solver"].GetString()]
        AddEchoLevelToSettings(predictor_settings, parent_echo_level)
        predictors.append(CreatePredictor(predictor_settings, solver))
    return predictors

def CreateConvergenceAccelerators(convergence_accelerator_settings_list, solvers, parent_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator
    convergence_accelerators = []
    for conv_acc_settings in convergence_accelerator_settings_list:
        solver = solvers[conv_acc_settings["solver"].GetString()]
        AddEchoLevelToSettings(conv_acc_settings, parent_echo_level)
        convergence_accelerators.append(CreateConvergenceAccelerator(conv_acc_settings, solver))

    return convergence_accelerators

def CreateConvergenceCriteria(convergence_criteria_settings_list, solvers, parent_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.convergence_criteria_factory import CreateConvergenceCriteria
    convergence_criteria = []
    for conv_crit_settings in convergence_criteria_settings_list:
        solver = solvers[conv_crit_settings["solver"].GetString()]
        AddEchoLevelToSettings(conv_crit_settings, parent_echo_level)
        convergence_criteria.append(CreateConvergenceCriteria(conv_crit_settings, solver))

    return convergence_criteria

def CreateCouplingOperations(coupling_operations_settings_dict, solvers, parent_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.coupling_operation_factory import CreateCouplingOperation
    coupling_operations = {}
    for coupling_operation_name, coupling_operation_settings in coupling_operations_settings_dict.items():
        AddEchoLevelToSettings(coupling_operation_settings, parent_echo_level)
        coupling_operations[coupling_operation_name] = CreateCouplingOperation(coupling_operation_settings, solvers)

    return coupling_operations

def CreateDataTransferOperators(data_transfer_operators_settings_dict, parent_echo_level):
    from KratosMultiphysics.CoSimulationApplication.factories.data_transfer_operator_factory import CreateDataTransferOperator
    data_transfer_operators = {}
    for data_transfer_operators_name, data_transfer_operators_settings in data_transfer_operators_settings_dict.items():
        AddEchoLevelToSettings(data_transfer_operators_settings, parent_echo_level)
        data_transfer_operators[data_transfer_operators_name] = CreateDataTransferOperator(data_transfer_operators_settings)

    return data_transfer_operators
