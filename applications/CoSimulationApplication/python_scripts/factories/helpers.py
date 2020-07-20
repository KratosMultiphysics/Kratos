# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.coupling_operation_factory import CreateCouplingOperation
from KratosMultiphysics.CoSimulationApplication.factories.data_transfer_operator_factory import CreateDataTransferOperator
from KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper import ConvergenceAcceleratorWrapper
from KratosMultiphysics.CoSimulationApplication.convergence_criteria.convergence_criteria_wrapper import ConvergenceCriteriaWrapper
from KratosMultiphysics.CoSimulationApplication.factories.predictor_factory import CreatePredictor


def AddEchoLevelToSettings(settings, echo_level):
    echo_level_params = KM.Parameters("""{
        "echo_level" : %d
    }""" % echo_level)
    settings.AddMissingParameters(echo_level_params)


def CreatePredictors(predictor_settings_list, solvers, parent_echo_level):
    predictors = []
    for predictor_settings in predictor_settings_list:
        solver = solvers[predictor_settings["solver"].GetString()]
        AddEchoLevelToSettings(predictor_settings, parent_echo_level)
        predictors.append(CreatePredictor(predictor_settings, solver))
    return predictors

def CreateConvergenceAccelerators(convergence_accelerator_settings_list, solvers, parent_echo_level):
    convergence_accelerators = []
    for conv_acc_settings in convergence_accelerator_settings_list:
        solver = solvers[conv_acc_settings["solver"].GetString()]
        AddEchoLevelToSettings(conv_acc_settings, parent_echo_level)
        convergence_accelerators.append(ConvergenceAcceleratorWrapper(conv_acc_settings, solver))

    return convergence_accelerators

def CreateConvergenceCriteria(convergence_criterion_settings_list, solvers, parent_echo_level):
    convergence_criteria = []
    for conv_crit_settings in convergence_criterion_settings_list:
        solver = solvers[conv_crit_settings["solver"].GetString()]
        AddEchoLevelToSettings(conv_crit_settings, parent_echo_level)
        convergence_criteria.append(ConvergenceCriteriaWrapper(conv_crit_settings, solver))

    return convergence_criteria

def CreateCouplingOperations(coupling_operations_settings_dict, solvers, parent_coupled_solver_process_info, parent_echo_level):
    coupling_operations = {}
    for coupling_operation_name, coupling_operation_settings in coupling_operations_settings_dict.items():
        AddEchoLevelToSettings(coupling_operation_settings, parent_echo_level)
        coupling_operations[coupling_operation_name] = CreateCouplingOperation(coupling_operation_settings, solvers, parent_coupled_solver_process_info)

    return coupling_operations

def CreateDataTransferOperators(data_transfer_operators_settings_dict, parent_echo_level):
    data_transfer_operators = {}
    for data_transfer_operators_name, data_transfer_operators_settings in data_transfer_operators_settings_dict.items():
        AddEchoLevelToSettings(data_transfer_operators_settings, parent_echo_level)
        data_transfer_operators[data_transfer_operators_name] = CreateDataTransferOperator(data_transfer_operators_settings)

    return data_transfer_operators
