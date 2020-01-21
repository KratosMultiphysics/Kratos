from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.coupling_operation_factory import CreateCouplingOperation
from KratosMultiphysics.CoSimulationApplication.factories.data_transfer_operator_factory import CreateDataTransferOperator
from KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper import ConvergenceAcceleratorWrapper
from KratosMultiphysics.CoSimulationApplication.convergence_criteria.convergence_criteria_wrapper import ConvergenceCriteriaWrapper
from KratosMultiphysics.CoSimulationApplication.factories.predictor_factory import CreatePredictor


import KratosMultiphysics.CoSimulationApplication.colors as colors

### This file contains functionalities that are commonly used in CoSimulation ###

def cs_print_info(label, *args):
    KM.Logger.PrintInfo(colors.bold(label), " ".join(map(str,args)))

def cs_print_warning(label, *args):
    KM.Logger.PrintWarning(colors.bold(label), " ".join(map(str,args)))


def UsingPyKratos():
    return any(["pyKratos" in i_path for i_path in KM.__path__])

def SettingsTypeCheck(settings):
    if not isinstance(settings, KM.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")


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

def CreateCouplingOperations(coupling_operations_settings_dict, solvers, parent_echo_level):
    coupling_operations = {}
    for coupling_operation_name, coupling_operation_settings in coupling_operations_settings_dict.items():
        AddEchoLevelToSettings(coupling_operation_settings, parent_echo_level)
        coupling_operations[coupling_operation_name] = CreateCouplingOperation(coupling_operation_settings, solvers)

    return coupling_operations

def CreateDataTransferOperators(data_transfer_operators_settings_dict, parent_echo_level):
    data_transfer_operators = {}
    for data_transfer_operators_name, data_transfer_operators_settings in data_transfer_operators_settings_dict.items():
        AddEchoLevelToSettings(data_transfer_operators_settings, parent_echo_level)
        data_transfer_operators[data_transfer_operators_name] = CreateDataTransferOperator(data_transfer_operators_settings)

    return data_transfer_operators

def AllocateHistoricalVariablesFromCouplingData(data_list, model, solver_name):
    '''This function retrieves the historical variables that are needed for the ModelParts from the
    specified CouplingInterfaceDatas and allocates them on the ModelParts
    Note that it can only be called after the (Main-)ModelParts are created
    '''
    for data in data_list:
        hist_var_dict = data.GetHistoricalVariableDict()
        for full_model_part_name, variable in hist_var_dict.items():
            main_model_part_name = full_model_part_name.split(".")[0]
            if not model.HasModelPart(main_model_part_name):
                raise Exception('ModelPart "{}" does not exist in solver "{}"!'.format(main_model_part_name, solver_name))
            main_model_part = model[main_model_part_name]
            if not main_model_part.HasNodalSolutionStepVariable(variable):
                cs_print_info("CoSimTools", 'Allocating historical variable "{}" in ModelPart "{}" for solver "{}"'.format(variable.Name(), main_model_part_name, solver_name))
                main_model_part.AddNodalSolutionStepVariable(variable)

def CreateMainModelPartsFromCouplingData(data_list, model, solver_name):
    '''This function creates the Main-ModelParts that are used in the specified CouplingInterfaceDatas
    '''
    for data in data_list:
        main_model_part_name = data.model_part_name.split(".")[0]
        if not model.HasModelPart(main_model_part_name):
            model.CreateModelPart(main_model_part_name)
            cs_print_info("CoSimTools", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, solver_name))

def RecursiveCreateModelParts(model_part, model_part_name):
    '''This function creates a hierarchy of SubModelParts on a given ModelPart
    '''
    model_part_name, *sub_model_part_names = model_part_name.split(".")
    if not model_part.HasSubModelPart(model_part_name):
        cs_print_info("CoSimTools", 'Created "{}" as SubModelPart of "{}"'.format(model_part_name, model_part.Name))
        model_part = model_part.CreateSubModelPart(model_part_name)
    if len(sub_model_part_names) > 0:
        RecursiveCreateModelParts(model_part, ".".join(sub_model_part_names))
