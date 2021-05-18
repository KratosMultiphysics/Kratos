# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import cs_print_info
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

def AllocateHistoricalVariablesFromCouplingDataSettings(data_settings_list, model, solver_name):
    '''This function allocates historical variables for the Modelparts

    It retrieves the historical variables that are needed for the ModelParts from the
    specified CouplingInterfaceData-settings and allocates them on the ModelParts
    Note that it can only be called after the (Main-)ModelParts are created
    '''
    data_settings_list = data_settings_list.Clone() # clone to not mess with the following validation

    for data_settings in data_settings_list.values():
        CouplingInterfaceData.GetDefaultParameters()
        data_settings.ValidateAndAssignDefaults(CouplingInterfaceData.GetDefaultParameters())

        if data_settings["location"].GetString() == "node_historical":
            variable = KM.KratosGlobals.GetVariable(data_settings["variable_name"].GetString())
            main_model_part_name = data_settings["model_part_name"].GetString().split(".")[0]
            if not model.HasModelPart(main_model_part_name):
                raise Exception('ModelPart "{}" does not exist in solver "{}"!'.format(main_model_part_name, solver_name))
            main_model_part = model[main_model_part_name]
            if not main_model_part.HasNodalSolutionStepVariable(variable):
                cs_print_info("CoSimTools", 'Allocating historical variable "{}" in ModelPart "{}" for solver "{}"'.format(variable.Name(), main_model_part_name, solver_name))
                main_model_part.AddNodalSolutionStepVariable(variable)

def CreateMainModelPartsFromCouplingDataSettings(data_settings_list, model, solver_name):
    '''This function creates the Main-ModelParts that are used in the specified CouplingInterfaceData-settings'''
    data_settings_list = data_settings_list.Clone() # clone to not mess with the following validation

    for data_settings in data_settings_list.values():
        CouplingInterfaceData.GetDefaultParameters()
        data_settings.ValidateAndAssignDefaults(CouplingInterfaceData.GetDefaultParameters())

        main_model_part_name = data_settings["model_part_name"].GetString().split(".")[0]
        if not model.HasModelPart(main_model_part_name):
            model.CreateModelPart(main_model_part_name)
            cs_print_info("CoSimTools", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, solver_name))

def RecursiveCreateModelParts(model_part, model_part_name):
    '''This function creates a hierarchy of SubModelParts on a given ModelPart'''
    model_part_name, *sub_model_part_names = model_part_name.split(".")
    if model_part.HasSubModelPart(model_part_name):
        model_part = model_part.GetSubModelPart(model_part_name)
    else:
        cs_print_info("CoSimTools", 'Created "{}" as SubModelPart of "{}"'.format(model_part_name, model_part.Name))
        model_part = model_part.CreateSubModelPart(model_part_name)
    if len(sub_model_part_names) > 0:
        RecursiveCreateModelParts(model_part, ".".join(sub_model_part_names))

def CreateModelPartsFromCouplingDataSettings(data_settings_list, model, solver_name):
    '''This function creates the ModelParts-hierarchie that are used in the specified CouplingInterfaceData-settings'''
    data_settings_list = data_settings_list.Clone() # clone to not mess with the following validation

    for data_settings in data_settings_list.values():
        CouplingInterfaceData.GetDefaultParameters()
        data_settings.ValidateAndAssignDefaults(CouplingInterfaceData.GetDefaultParameters())

        splitted_name = data_settings["model_part_name"].GetString().split(".")
        main_model_part_name = splitted_name[0]
        sub_model_part_names = splitted_name[1:]
        if model.HasModelPart(main_model_part_name):
            main_model_part = model.GetModelPart(main_model_part_name)
        else:
            main_model_part = model.CreateModelPart(main_model_part_name)
            cs_print_info("CoSimTools", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, solver_name))
        if len(sub_model_part_names) > 0:
            RecursiveCreateModelParts(main_model_part, ".".join(sub_model_part_names))
